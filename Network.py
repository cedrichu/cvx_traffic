import numpy as np
from mpi4py import MPI
from cvxpy import *
import Solver 
import Constants
from scipy.optimize import fsolve

class TrafficNetwork(object):
	def __init__(self, agent_list, adjacent_matrix):
		self._vertex = agent_list
		self._vertex_num = len(agent_list)
		self._adjacent_matrix = adjacent_matrix
		
		self.connect_queues()
		self.set_turn_prop()
		self.set_turn_prop_up()

	def connect_queues(self):
		for vi in range(self._vertex_num):
			for vj in range(self._vertex_num):
				if self._adjacent_matrix[vi][vj] != []:
					self._vertex[vi].set_neighbors(self._vertex[vj], self._adjacent_matrix[vi][vj])

	def set_turn_prop(self):
		for v in self._vertex:
			v.set_turn_prop()

	def set_turn_prop_up(self):
		for v in self._vertex:
			v.set_turn_prop_up()


	def __repr__(self):
		return str([v._agent_id for v in self._vertex])


class TrafficAgentModel(object):
	def __init__(self, agent_id, local_queue_num):
		self._agent_id = agent_id
		self._local_queue_num = local_queue_num
		self._local_queue = []
		self._neighbor_agents = []

		self.init_local_queue()

		'''constants'''
		self._sat_flow_rate = Constants.SAT_FLOW_RATE

	def init_local_queue(self):
		for i in range(self._local_queue_num):
			queue = TrafficQueue(self._agent_id, i)
			self._local_queue.append(queue)

	def init_queue_solver_vars(self):		
		for i in range(self._local_queue_num):
			self.get_local_queue(i).Initialize_variables()			

	def set_neighbors(self, neighbor, edge):
		
		if neighbor not in self._neighbor_agents:
			self._neighbor_agents.append(neighbor)

		for q in edge[0]:
			self._local_queue[q].set_downstream(neighbor.get_local_queue(edge[1]))
			neighbor.get_local_queue(edge[1]).set_upstream(self._local_queue[q])

	def set_ext_arr_rate(self, rate):
		for idx,q  in enumerate(self._local_queue):
			q.set_ext_arr_rate(rate[idx])

	def set_turn_prop(self):
		for v in self._local_queue:
			v.set_turn_prop()

	def set_turn_prop_up(self):
		for v in self._local_queue:
			v.set_turn_prop_up()

	def get_local_queue(self, queue_id):
		return self._local_queue[queue_id]

	def get_agent_id(self):
		return self._agent_id

	def __repr__(self):
		return str([self._agent_id, [v.get_agent_id() for v in self._neighbor_agents]])

	def get_queue_sum_sat_constraint(self):
		eqn = Variable()
		
		for i in range(self._local_queue_num):
			eqn = eqn + self.get_local_queue(i)._vars[4][0]
		
		return [eqn == self._sat_flow_rate]
		
	def get_all_constraints(self):
		constraints = []

		for i in range(self._local_queue_num):
			constraints += self.get_local_queue(i).get_constraints()	

		constraints += self.get_queue_sum_sat_constraint()	
		return constraints	

	def get_new_objective(self):		
		obj = 0
		
		for i in range(self._local_queue_num):
			obj = obj + self.get_local_queue(i).get_objective()			
		
		return obj	

	def get_primal_objective(self):
		obj = 0
		for i in range(self._local_queue_num):
			obj = obj + self.get_local_queue(i).get_primal_obj()
		return obj

	def compute_residuals(self):	
		obj = 0
		for i in range(self._local_queue_num):
			obj = obj + self.get_local_queue(i).compute_residuals()
		return obj

	def get_new_data(self):
		comm = MPI.COMM_WORLD
		lb = comm.recv(source = Constants.BB_TREE_ID , tag = Constants.NEW_LOWER_BOUNDS)	
		ub = comm.recv(source = Constants.BB_TREE_ID , tag = Constants.NEW_UPPER_BOUNDS)
		dual = comm.recv(source = Constants.BB_TREE_ID , tag = Constants.NEW_DUAL_VARS)
		self.set_new_data(lb, ub, dual)

	def set_new_data(self, lb , ub , dual):	
		self.lb = lb
		self.ub = ub
		self.dual = dual

	def construct_prob(self):
		comm = MPI.COMM_WORLD
		comm.send(lb , dest = Constants.BB_TREE_ID, tag = Constants.TREE_START_LB)	
		comm.send(ub , dest = Constants.BB_TREE_ID, tag = Constants.TREE_START_UB)
		comm.send(dual , dest = Constants.BB_TREE_ID, tag = Constants.TREE_START_DUAL)

	def Update_consensus_vars(self):
		for i in range(self._local_queue_num):
			self.get_local_queue(i).send_rel_vars()

		neigh_data = dict()

		for i in range(self._local_queue_num):
			neigh_data[i] = self.get_local_queue(i).receive_rel_vars()

		for i in range(self._local_queue_num):	
			self.get_local_queue(i).solve_send_vars(neigh_data[i])

		for i in range(self._local_queue_num):	
			self.get_local_queue(i).recv_update_solved_vars()			

	def get_consensus_constraints(self):
		constraints = []
		for i in range(self._local_queue_num):
			constraints += self.get_local_queue(i).get_consensus_constraints()
		return constraints	

	def get_coupling_constraints(self , agent_list):
		constraints = []
		for i in range(self._local_queue_num):
			constraints += self.get_local_queue(i).get_coupling_constraints(agent_list)
		return constraints		

	def Update_Dual_Vars(self):
		for i in range(self._local_queue_num):
			self.get_local_queue(i).Update_Dual_Vars()

	def Check_Feasibility(self):
		for i in range(self._local_queue_num):	
			self.get_local_queue(i).feasibility()

	def solve_problems(self):	 	
	 	constraints =  self.get_all_constraints()	 	
	 	
	 	comm = MPI.COMM_WORLD
		rank = comm.Get_rank()
	 	
	 	iter = 0

	 	while(iter < 200):	
	 		obj = self.get_new_objective()	
			prob = Problem(Minimize(obj), constraints)
			prob.solve(solver=ECOS , warm_start = True , max_iters = 2000 , abstol = 10 ** -11)
			#print prob.value
			self.Update_consensus_vars()
			
			#print self.compute_residuals() , iter
			self.Update_Dual_Vars()	 	
			
			print 'feasibility agent_id', self._agent_id
			print self.get_primal_objective().value , iter
			
			self.Check_Feasibility()


			iter = iter + 1
	
class TrafficQueue(object):
	def __init__(self, agent_id, queue_id):
		self._agent_id = agent_id
		self._queue_id = queue_id
		self.lb = []
		self.ub = []
		self._constraints = []
		self._vars = []
		self._dual_vars	= []
		self._coup_res = dict()			

		'''local variables'''
		self._arr_rate = 0
		self._eff_arr_rate = 0
		self._block_prob = 0
		self._service_rate = 0
		self._eff_service_rate = 0
		self._unblocking_rate = 0
		self_being_blocked_prob = 0
		self._ro = 0
		self._v = 0
		self._d = 0
		self._t = 0
		self._w = 0
		self._g = 0
		self._omega = 0
		self._bz = 0
		self._cz = 0
		self._ez = 0
		
		'''neighbors'''
		self._upstream_queue = []
		self._downstream_queue = []
		self._az = []
		self._dz = []
		self._fz = []


		'''constants'''
		self._turn_prop = {}
		self._turn_prop_up = {}
		self._ext_arr_rate = 0
		self._capacity = Constants.CAPACITY
		self._epsilon = Constants.EPSILON
		self._speed_limit = Constants.SPEED_LIMIT

	def init_lb(self):
		self.lb.append([0]) #lambda
		self.lb.append([self._epsilon]) #lambda^eff
		self.lb.append([0]) #P(N = k)
		self.lb.append([self._epsilon]) #mu_eff
		self.lb.append([self._epsilon]) #mu
		self.lb.append([self._epsilon]) #mu^tilde
		self.lb.append([0]) #P
		self.lb.append([self._epsilon]) #rho
		self.lb.append([0]) #v
		
		self.lb.append([0]) #d
		self.lb.append([0]) #t
		self.lb.append([self._epsilon]) #w
		self.lb.append([0]) #f
		self.lb.append([0]) #r


		lb_consus_down1= {}
		for v in self._downstream_queue:
			lb_consus_down1[(v.get_agent_id(), v.get_queue_id())] = self._epsilon 
		self.lb.append(lb_consus_down1) #A
		self.lb.append([-self._speed_limit]) #B
		self.lb.append([self._epsilon]) #C
		lb_consus_down2= {}
		for v in self._upstream_queue:
			lb_consus_down2[(v.get_agent_id(), v.get_queue_id())] = self.lb[12][0] 
		self.lb.append(lb_consus_down2) #D
		self.lb.append([0]) #E
		lb_consus_down3= {}
		for v in self._upstream_queue:
			lb_consus_down3[(v.get_agent_id(), v.get_queue_id())] = 0
		self.lb.append(lb_consus_down3) #F
		self.lb.append([1]) #d_i\mu_i^eff
		self.lb.append([1]) #w_i\mu_i

	def init_ub(self):
		self.ub.append([self._speed_limit])#lambda
		self.ub.append([self._speed_limit])#lambda^eff
		self.ub.append([1])#P(N = k)
		self.ub.append([self._speed_limit])#mu_eff
		self.ub.append([self._speed_limit])#mu
		self.ub.append([1])#mu^tilde
		self.ub.append([1])#P
		self.ub.append([self._speed_limit/self._epsilon])#rho
		self.ub.append([1])#v



		self.ub.append([2/self._epsilon])#d
		self.ub.append([2/self._epsilon])#t
		self.ub.append([2/self._epsilon])#w
		self.ub.append([2*self._speed_limit/self._epsilon])#f
		self.ub.append([2*self._speed_limit/self._epsilon])#r


		ub_consus_down1= {}
		for v in self._downstream_queue:
			ub_consus_down1[(v.get_agent_id(), v.get_queue_id())] = self._speed_limit*self._turn_prop[(v.get_agent_id(), v.get_queue_id())]
		self.ub.append(ub_consus_down1) #A
		self.ub.append([self._speed_limit]) #B
		self.ub.append([2*len(self._downstream_queue)*self._speed_limit/self._epsilon]) #C
		ub_consus_down2= {}
		for v in self._upstream_queue:
			ub_consus_down2[(v.get_agent_id(), v.get_queue_id())] = self._speed_limit/self._epsilon
		self.ub.append(ub_consus_down2) #D
		self.ub.append([1]) #E
		ub_consus_down3= {}
		for v in self._upstream_queue:
			ub_consus_down3[(v.get_agent_id(), v.get_queue_id())] = self._turn_prop_up[(v.get_agent_id(), v.get_queue_id())]
		self.ub.append(ub_consus_down3) #F
		self.ub.append([1]) #d_i\mu_i^eff
		self.ub.append([1]) #w_i\mu_i

	def Initialize_variables(self):
		self._vars , self._dual_vars = Solver.Variable_Initialization(self._upstream_queue , self._downstream_queue)
		self.init_lb()
		self.init_ub()	

	def set_upstream(self, queue):
		self._upstream_queue.append(queue)

	def set_downstream(self, queue):
		self._downstream_queue.append(queue)

	def set_turn_prop(self, prop=[0.1,0.8,0.1]):
		for idx, v in enumerate(self._downstream_queue):
			self._turn_prop[(v._agent_id, v._queue_id)] = prop[idx]

	def set_turn_prop_up(self):
		for idx, v in enumerate(self._upstream_queue):
			self._turn_prop_up[(v._agent_id, v._queue_id)] = v.get_turn_prop(self._agent_id, self._queue_id)

	def set_ext_arr_rate(self, rate):
		self._ext_arr_rate = rate

	def set_capacity(self, capacity):
		self._capacity = capacity

	def get_agent_id(self):
		return self._agent_id

	def get_queue_id(self):
		return self._queue_id

	def get_turn_prop(self, agent_id, queue_id):
		return self._turn_prop[(agent_id, queue_id)]

	def get_turn_prop_up(self, agent_id, queue_id):
		return self._turn_prop_up[(agent_id, queue_id)]

	def __repr__(self):
		return str(self._agent_id)+' '+str(self._queue_id)

	def create_constraints(self):
	 	return Solver.Create_constraints_for_queue(self._vars , self._upstream_queue, self._downstream_queue, self.lb, self.ub , self._capacity)

	def get_constraints(self):
		_constraints = self.create_constraints()
	 	return _constraints	

	def get_objective(self):
		return Solver.Get_total_queue_objective(self._vars, self._dual_vars , self._upstream_queue , self._downstream_queue , self._ext_arr_rate, self._turn_prop , self._turn_prop_up) 	
		
	def Update_Dual_Vars(self):
		return Solver.Update_Dual_Vars(self._vars, self._dual_vars , self._upstream_queue , self._downstream_queue , self._ext_arr_rate, self._turn_prop, self._turn_prop_up , self._coup_res)

	def send_rel_vars(self):
		return Solver.send_rel_vars(self._vars , self._dual_vars , self._upstream_queue, self._downstream_queue, self._turn_prop , self._queue_id , self.lb , self.ub , self._turn_prop_up)
		
	def receive_rel_vars(self):
		return Solver.receive_rel_vars(self._vars , self._dual_vars , self._upstream_queue , self._downstream_queue , self._queue_id )

	def solve_send_vars(self , neigh_data):
		return Solver.solve_coupling_eqns_send_sols(self._vars , self._dual_vars , self._upstream_queue , self._downstream_queue , self._queue_id , neigh_data , self._ext_arr_rate , self.lb , self.ub , self._coup_res)

	def recv_update_solved_vars(self):
		return Solver.recv_update_solved_vars(self._vars , self._agent_id , self._queue_id , self._upstream_queue , self._downstream_queue)

	def get_primal_obj(self):
		obj = inv_pos(1 - self._vars[7][0]) - 1
		return obj

	def get_consensus_constraints(self):
		return Solver.get_consensus_constraints(self._vars, self._dual_vars , self._upstream_queue , self._downstream_queue , self._ext_arr_rate, self._turn_prop , self._turn_prop_up)

	def compute_residuals(self):
		return Solver.compute_residuals(self._vars, self._dual_vars , self._upstream_queue , self._downstream_queue , self._ext_arr_rate, self._turn_prop , self._turn_prop_up , self._coup_res)	

	def get_coupling_constraints(self , agent_list):
		return Solver.get_coupling_constraints(self._vars, self._upstream_queue , self._downstream_queue , agent_list , self._agent_id , self._queue_id)

	def feasibility(self):

		zG = []
		for i in range(14):
			zG.append(self._vars[i][0].value)


		solution,infodict,ier, msg = fsolve(self.system_equations, zG, full_output=True,maxfev=100,xtol=0.1)
		
		print 'feasibility'
		print solution
		#print infodict
		print ier
		print msg
		
		return ier

	
	def system_equations(self, var):
		F = []
		
		F += [var[0] - var[1]*var[8]] #(1)
		
		x = 0
		for q in self._downstream_queue: 
			x += pow(self._vars[14][(q._agent_id, q._queue_id)].value  - self._turn_prop[(q._agent_id,q._queue_id)]*var[1],2) #(2)
		x += pow(self._vars[15].value - var[1] + self._ext_arr_rate*(1-var[2]),2) #(2)
		F += [x]
		
		F += [var[9] - var[11] - var[6]*var[10]]  #(3)
		

		y = 0
		y += pow(self._vars[16].value - var[1]*var[10],2) #(4)
		for q in self._upstream_queue:
			y += pow(self._vars[17][(q._agent_id, q._queue_id)].value - var[1]*var[9],2)#(4)
		F += [y]
		
		z = 0
		z += pow(self._vars[18].value - var[6],2)#(5)
		for q in self._upstream_queue:
			z += pow(self._vars[19][(q._agent_id, q._queue_id)].value - self._turn_prop_up[(q._agent_id,q._queue_id)]*var[2],2) #(5)
		F += [z]
		
		k = pow(var[7],self._capacity)
		F += [var[2] - (k-pow(var[7],self._capacity+1))]#/(1-k*var[7])] #6

		F += [var[7] - var[0]*var[9]] #7

		F += [var[8] - 1 + var[2]]
		F += [var[1] - var[0]*var[8]] #8
		F += [var[9]*var[3]  - 1] #9
		F += [var[11]*var[4]  - 1] #10
		F += [var[10]*var[5]  - var[6]] #11
		F += [var[1] - var[13]*var[5]] #12
		F += [var[1] - var[12] * var[3]] #13
		return F

	def system_equations2(self):
		F = []
		
		F += [var[0] - var[1]*var[8]] #(1)
		
		x = 0
		# for q in self._downstream_queue: 
		# 	x += pow(self._vars[14][(q._agent_id, q._queue_id)].value  - self._turn_prop[(q._agent_id,q._queue_id)]*var[1],2) #(2)
		x += self._vars[15].value - var[1] + self._ext_arr_rate*(1-var[2]) #(2)
		F += [x]
		
		F += [var[9] - var[11] - var[6]*var[10]]  #(3)
		

		y = 0
		y += self._vars[16].value - var[1]*var[10] #(4)
		# for q in self._upstream_queue:
		# 	y += pow(self._vars[17][(q._agent_id, q._queue_id)].value - var[1]*var[9],2)#(4)
		F += [y]
		
		z = 0
		z += self._vars[18].value - var[6]#(5)
		# for q in self._upstream_queue:
		# 	z += pow(self._vars[19][(q._agent_id, q._queue_id)].value - self._turn_prop_up[(q._agent_id,q._queue_id)]*var[2],2) #(5)
		F += [z]
		
		k = pow(var[7],self._capacity)
		F += [var[2] - (k-pow(var[7],self._capacity+1))]#/(1-k*var[7])] #6

		F += [var[7] - var[0]*var[9]] #7

		F += [var[8] - 1 + var[2]]
		F += [var[1] - var[0]*var[8]] #8
		F += [var[9]*var[3]  - 1] #9
		F += [var[11]*var[4]  - 1] #10
		F += [var[10]*var[5]  - var[6]] #11
		F += [var[1] - var[13]*var[5]] #12
		F += [var[1] - var[12] * var[3]] #13
		return F


# def main():
# 	'''create agents'''
# 	queue_num = 4
# 	t1 = TrafficAgentModel(0,queue_num)
# 	t2 = TrafficAgentModel(1,queue_num)
# 	t3 = TrafficAgentModel(2,queue_num)
# 	t4 = TrafficAgentModel(3,queue_num)
# 	agent_list = [t1,t2,t3,t4]
# 	'''specify external arrival rates'''
# 	t1.set_ext_arr_rate([10,0,0,10])
# 	t2.set_ext_arr_rate([10,10,0,0])
# 	t3.set_ext_arr_rate([0,10,10,0])
# 	t4.set_ext_arr_rate([0,0,10,10])
# 	'''specify connection'''
# 	n = len(agent_list)
# 	adjacent_matrix = [[[]for x in range(n)] for y in range(n)] 
# 	adjacent_matrix[0][1] = [[0,3,2],3]
# 	adjacent_matrix[0][3] = [[1,0,3],0]
# 	adjacent_matrix[3][2] = [[0,3,2],3]
# 	adjacent_matrix[3][0] = [[1,2,3],2]
# 	adjacent_matrix[1][2] = [[1,0,3],0]
# 	adjacent_matrix[1][0] = [[0,1,2],1]
# 	adjacent_matrix[2][1] = [[1,2,3],2]
# 	adjacent_matrix[2][3] = [[0,1,2],1]
# 	'''construct network'''
# 	network = TrafficNetwork(agent_list, adjacent_matrix)

# 	'''for testing'''
# 	print adjacent_matrix
# 	print network
# 	for t in agent_list:
# 		print t
# 	for q in t1._local_queue:
# 		print q
	



if __name__ == '__main__':
	main()




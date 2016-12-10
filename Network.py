import numpy as np
from mpi4py import MPI
from cvxpy import *
 

class TrafficNetwork(object):
	def __init__(self, agent_list, adjacent_matrix):
		self._vertex = agent_list
		self._vertex_num = len(agent_list)
		self._adjacent_matrix = adjacent_matrix
		
		self.connect_queues()
		self.set_turn_prop()

	def connect_queues(self):
		for vi in range(self._vertex_num):
			for vj in range(self._vertex_num):
				if self._adjacent_matrix[vi][vj] != []:
					self._vertex[vi].set_neighbors(self._vertex[vj], self._adjacent_matrix[vi][vj])

	def set_turn_prop(self):
		for v in self._vertex:
			v.set_turn_prop()

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
		self._sat_flow_rate = 0

	def init_local_queue(self):
		for i in range(self._local_queue_num):
			queue = TrfficQueue(self._agent_id, i)
			self._local_queue.append(queue)

	def set_neighbors(self, neighbor, edge):
		
		if neighbor not in self._neighbor_agents:
			self._neighbor_agents.append(neighbor)

		for q in edge[0]:
			self._local_queue[edge[1]].set_upstream(neighbor.get_local_queue(q))
			neighbor.get_local_queue(q).set_downstream(self._local_queue[edge[1]])

	def set_ext_arr_rate(self, rate):
		for idx,q  in enumerate(self._local_queue):
			q.set_ext_arr_rate(rate[idx])

	def set_turn_prop(self):
		for v in self._local_queue:
			v.set_turn_prop()

	def get_local_queue(self, queue_id):
		return self._local_queue[queue_id]

	def get_agent_id(self):
		return self._agent_id

	def __repr__(self):
		return str([self._agent_id, [v.get_agent_id() for v in self._neighbor_agents]])

	def get_queue_sum_sat_constraint(self):
		eqn = Variable()
		constraint = []

		for i in range(self._local_queue_num)
			eqn = eqn + self.get_local_queue(i)._vars[4][0]
		
		constraint.append(eqn == self._sat_flow_rate)
		return constraint

	def get_all_constraints(self):
		constraints = []

		for i in range(self._local_queue_num)
			constraints.append(self.get_local_queue(i).get_constraints())	

		constraints.append(self.get_queue_sum_sat_constraint())	
		return constraints	

	def get_new_objective(self):		
		obj = Minimize(pos(0))	

		for i in range(self._local_queue_num)
			obj = obj + self.get_local_queue(i).get_objective()

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

	#updates only dual variables corresponding to consensus constraints, 
	#dual variables w.r.t coupling constraints are updated in previous solve_send_vars() of Update_consensus_vars() 
	def Update_Dual_Vars(self):
		for i in range(self._local_queue_num):
			self.get_local_queue(i).Update_Dual_Vars()

	def solve_problems(self):	 	
	 	constraints =  self.get_all_constraints()	 	
	 	while(1):	
	 		obj = self.get_new_objective()
	 		prob = Problem(Minimize(obj), constraints)
			prob.solve()

			self.Update_consensus_vars()
			self.Update_Dual_Vars()	 		
	
class TrfficQueue(object):
	def __init__(self, agent_id, queue_id):
		self._agent_id = agent_id
		self._queue_id = queue_id

		self._vars , self._dual_vars = Solver.Variable_Initialization()
		self._constraints = []
		self.lb = []
		self.ub = []
						
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
		self._ext_arr_rate = 0
		self._capacity = 20
		self._epsilon = 0.01
		self._speed_limit = 20.0

	def init_lb(self):
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])

		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])
		self.lb.append([self._epsilon])


		lb_consus_down1= {}
		for v in self._downstream_queue:
			lb_consus_down1[(v.get_agent_id(), v.get_queue_id())] = self._epsilon 
		self.lb.append(lb_consus_down1)
		self.lb.append([-self._speed_limit])
		self.lb.append([self._epsilon])
		lb_consus_down2= {}
		for v in self._uptream_queue:
			lb_consus_down2[(v.get_agent_id(), v.get_queue_id())] = self._epsilon 
		self.lb.append(lb_consus_down2)
		self.lb.append([self._epsilon])
		lb_consus_down3= {}
		for v in self._uptream_queue:
			lb_consus_down3[(v.get_agent_id(), v.get_queue_id())] = self._epsilon 
		self.lb.append(lb_consus_down3)

	def init_ub(self):
		self.ub.append([self._speed_limit])
		self.ub.append([self._speed_limit])
		self.ub.append([1])
		self.ub.append([self._speed_limit])
		self.ub.append([self._speed_limit])
		self.ub.append([1])
		self.ub.append([1])
		self.ub.append([self._speed_limit/self._epsilon])

		self.ub.append([1/self._epsilon])
		self.ub.append([1/self._epsilon])
		self.ub.append([1/self._epsilon])
		self.ub.append([self._speed_limit/self._epsilon])
		self.ub.append([self._speed_limit/self._epsilon])


		ub_consus_down1= {}
		for v in self._downstream_queue:
			ub_consus_down1[(v.get_agent_id(), v.get_queue_id())] = self._speed_limit*self._turn_prop[(v.get_agent_id(), v.get_queue_id())]
		self.ub.append(ub_consus_down1)
		self.ub.append([self._speed_limit])
		self.ub.append([self._speed_limit/self._epsilon])
		ub_consus_down2= {}
		for v in self._upstream_queue:
			ub_consus_down2[(v.get_agent_id(), v.get_queue_id())] = self._speed_limit/self._epsilon
		self.ub.append(ub_consus_down2)
		self.ub.append([1])
		self.ub.append([1]*len(self._upstream_queue))
		ub_consus_down3= {}
		for v in self._upstream_queue:
			ub_consus_down3[(v.get_agent_id(), v.get_queue_id())] = 1
		self.ub.append(ub_consus_down3)



	def set_upstream(self, queue):
		self._upstream_queue.append(queue)

	def set_downstream(self, queue):
		self._downstream_queue.append(queue)

	def set_turn_prop(self, prop=[0.1,0.8,0.1]):
		for idx, v in enumerate(self._downstream_queue):
			self._turn_prop[(v._agent_id, v._queue_id)] = prop[idx]

	def set_ext_arr_rate(self, rate):
		self._ext_arr_rate = rate

	def set_capacity(self, capacity):
		self._capacity = capacity

	def get_agent_id(self):
		return self._agent_id

	def get_queue_id(self):
		return self._queue_id

	def get_turning_prop(self, agent_id, queue_id):
		return self._turn_prop[(agent_id, queue_id)]

	def __repr__(self):
		return str(self._agent_id)+' '+str(self._queue_id)

	def create_constraints(self):
	 	return Solver.Create_constraints_for_queue(self._vars , self._upstream_queue, self._downstream_queue, self.lb, self.ub , self._capacity)

	def get_constraints(self):
		_constraints = self.create_constraints()
	 	return _constraints	

	def get_objective(self):
		return Solver.Get_total_queue_objective(self._vars, self._dual_vars , self._upstream_queue , self._downstream_queue , self._ext_arr_rate, self._turn_prop) 	

	def Update_Dual_Vars(self):
		return Solver.Update_Dual_Vars(self._vars, self._dual_vars , self._upstream_queue , self._downstream_queue , self._ext_arr_rate, self._turn_prop)

	def send_rel_vars():
		return Solver.send_rel_vars(self._vars , self._dual_vars , self._upstream_queue, self._downstream_queue, self._turn_prop , self._queue_id )
		
	def receive_rel_vars():
		return Solver.receive_rel_vars(self._vars , self._dual_vars , self._upstream_queue , self._downstream_queue , self._queue_id )

	def solve_send_vars(neigh_data):
		return Solver.solve_coupling_eqns_send_sols(self._vars , self._dual_vars , self._upstream_queue , self._downstream_queue , self._queue_id , neigh_data , self._ext_arr_rate)

	def recv_update_solved_vars():
		return Solver.recv_update_solved_vars(self._vars , self._agent_id , self._queue_id , self._upstream_queue , self._downstream_queue)

def main():
	'''create agents'''
	queue_num = 4
	t1 = TrafficAgentModel(0,queue_num)
	t2 = TrafficAgentModel(1,queue_num)
	t3 = TrafficAgentModel(2,queue_num)
	t4 = TrafficAgentModel(3,queue_num)
	agent_list = [t1,t2,t3,t4]
	'''specify external arrival rates'''
	t1.set_ext_arr_rate([10,0,0,10])
	t2.set_ext_arr_rate([10,10,0,0])
	t3.set_ext_arr_rate([0,10,10,0])
	t4.set_ext_arr_rate([0,0,10,10])
	'''specify connection'''
	n = len(agent_list)
	adjacent_matrix = [[[]for x in range(n)] for y in range(n)] 
	adjacent_matrix[0][1] = [[0,3,2],3]
	adjacent_matrix[0][3] = [[1,0,3],0]
	adjacent_matrix[3][2] = [[0,3,2],3]
	adjacent_matrix[3][0] = [[1,2,3],2]
	adjacent_matrix[1][2] = [[1,0,3],0]
	adjacent_matrix[1][0] = [[0,1,2],1]
	adjacent_matrix[2][1] = [[1,2,3],2]
	adjacent_matrix[2][3] = [[0,1,2],1]
	'''construct network'''
	network = TrafficNetwork(agent_list, adjacent_matrix)

	'''for testing'''
	print adjacent_matrix
	print network
	for t in agent_list:
		print t
	for q in t1._local_queue:
		print q
	



if __name__ == '__main__':
	main()




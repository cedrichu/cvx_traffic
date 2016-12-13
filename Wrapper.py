from sortedcontainers import SortedDict
from mpi4py import MPI
from cvxpy import *
import math
import Constants
import Solver
import numpy as np
import Network
import matplotlib.pyplot as plt

class Traffic_Agent(object):
	#lb and ub should be a 2D list
	def __init__(self, neigh, iNumVars):
		self.upstream =  []
		self.downstream = []

		self.neigh = neigh	
		self.iNumVars = iNumVars
		self.Vars = Solver.Variable_Initialization(self, iNumVars)
		#Traffic_Agent.sample_constructor = Solver.Variable_Initialization

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

	def admm_optimize(self):
		#HCH - Fill
		return	
		
	def solve_problems(self):
		while(1):	
			self.get_new_data()	
			self.admm_optimize()

		

class Node(object):
	
	def __init__(self, lb, ub , dual):
		self.lb = lb
		self.ub = ub
		self.dual = dual
			

class BB_Tree_Agent(object):
	
	def __init__(self , iNumAgents):
		self.iNumAgents = iNumAgents
		#self.construct_prob()
		self.UB = float('inf')

	def construct_prob(self):
		comm = MPI.COMM_WORLD
		lb = dict()
		ub = dict()
		dual = dict()

		for i in range(0 , self.iNumAgents):
			lb[i] = comm.recv(source = i, tag = Constants.TREE_START_LB)
			ub[i] = comm.recv(source = i, tag = Constants.TREE_START_UB)
			dual[i] = comm.recv(source = i, tag = Constants.TREE_START_DUAL)				
		
		self.tree = SortedDict( [ (-float('inf') , Node( lb, ub , dual) ) ] )
			
	def get_next_node(self):
		return self.tree.popitem(last = False)

	def send_agent_node_info(self , node):
		comm = MPI.COMM_WORLD
		for i in range(0 , self.iNumAgents):
			data = node.lb[i]
			comm.send(data , dest = i , tag = Constants.NEW_LOWER_BOUNDS)
			data = node.ub[i]
			comm.send(data , dest = i , tag = Constants.NEW_UPPER_BOUNDS)
			data = node.dual[i]
			comm.send(data , dest = i , tag = Constants.NEW_DUAL_VARS)			
	
	def get_solution_from_agents(self):
		obj = 0
		feasible = True
		primal_sol = dict()
		dual_sol = dict()
		comm = MPI.COMM_WORLD
		for i in range(0 , self.iNumAgents):
			data = comm.recv(source = i , tag = Constants.GET_PRIMAL_SOLUTION)
			obj = obj + data[0];
			feasible = feasible * data[1]
			primal_sol[i] = data[2:]
			data = comm.recv(source = i , tag = Constants.GET_DUAL_VARS)
			dual_sol[i] = data
		return [obj , feasible, primal_sol, dual_sol]

	def create_new_nodes(self, obj, primal_sol , dual_sol , node):
		
		max = 0

		for i in range(0 , self.iNumAgents):	
			for j in range(0 , len(primal_sol[i])):
				for k in range(0 , len(primal_sol[i][j])):
					lbval = node.lb[i][j][k]
					ubval = node.ub[i][j][k]
					primal_val = primal_sol[i][j][k]
					ratio = ( (primal_val - lbval)/(ubval - lbval) ) * ( (ubval - primal_val)/(ubval - lbval) )

					if(ratio > max):
						split = [i , j , k]
						max = ratio

		if( 0 == max ):
			print 'unexpected behaviour in splitting'

		ub1 = node.ub
		ub1[split[0]][split[1]][split[2]] = primal_sol[i][j][k]	
			
		lb2 = node.lb
		lb2[split[0]][split[1]][split[2]] = primal_sol[i][j][k]

		self.tree[obj + random.random() * Constants.OBJ_PERT_FACTOR] = Node(node.lb, ub1 , dual_sol)
		self.tree[obj + random.random() * Constants.OBJ_PERT_FACTOR] = Node(lb2, ub , dual_sol)	
		
	def begin_optimization(self):
		iter = 0

		comm = MPI.COMM_WORLD
		obj_list = []
		inf_list = []
		x = []

		# datafile = open('data.txt','r')
		# for line in datafile:
		# 	data = line.split()
		# 	inf_list.append(float(data[1]))
		# 	x.append(int(data[2]))

		while(iter < Constants.MAX_ITER):
			obj = 0
			infeas = 0
			for i in range(0 , self.iNumAgents):
				data = comm.recv(source = i, tag = Constants.OPT_VAL)
				obj = obj + data[0]
				infeas = infeas + data[1] 

			print obj , math.sqrt(infeas) , iter	
			obj_list.append(obj)
			inf_list.append(math.sqrt(infeas))	
			x.append(iter)
			iter = iter + 1

		# plt.plot(np.array(x), np.array(obj_list))
		# #plt.legend( loc='upper right', numpoints = 1 )
		# plt.grid()
		# plt.ylabel('objective values')
		# plt.xlabel('iterations')
		# plt.show()

		# plt.plot(np.array(x), np.array(inf_list))
		# #plt.legend( loc='upper right', numpoints = 1 )
		# plt.grid()
		# plt.ylabel('residual')
		# plt.xlabel('iterations')
		# plt.show()

		

if __name__ == '__main__':
	
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	iNumAgents = 4;

	if ( Constants.BB_TREE_ID == rank ):
		tree = BB_Tree_Agent(iNumAgents)
		#node = tree.get_next_node()
		tree.begin_optimization()
	else:
		queue_num = 4
		t1 = Network.TrafficAgentModel(0,queue_num)
		t2 = Network.TrafficAgentModel(1,queue_num)
		t3 = Network.TrafficAgentModel(2,queue_num)
		t4 = Network.TrafficAgentModel(3,queue_num)
		agent_list = [t1,t2,t3,t4]
		'''specify external arrival rates'''
		ext_arr_rate = 5.0
		t1.set_ext_arr_rate([ext_arr_rate,0,0,ext_arr_rate])
		t2.set_ext_arr_rate([ext_arr_rate,ext_arr_rate,0,0])
		t3.set_ext_arr_rate([0,ext_arr_rate,ext_arr_rate,0])
		t4.set_ext_arr_rate([0,0,ext_arr_rate,ext_arr_rate])
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
		network = Network.TrafficNetwork(agent_list, adjacent_matrix)

		t = 0

		if(0 ==  rank):
			t = t1
		elif(1 == rank):
			t = t2
		elif(2 == rank):
			t = t3
		else:		
			t = t4

		t.init_queue_solver_vars()	
		t.solve_problems()




		




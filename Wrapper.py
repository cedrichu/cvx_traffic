from sortedcontainers import SortedDict
from mpi4py import MPI
from cvxpy import *
import math
import Constants
import Solver
import numpy

class Traffic_Agent(object):
	#lb and ub should be a 2D list
	def __init__(self, neigh, iNumVars):
		#self.upstream = 
		#self.downstream = 

		self.neigh = neigh	
		self.iNumVars = iNumVars
		self.Vars = []
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
		self.construct_prob()
		self.UB = float('inf')

	def construct_prob(self):
		comm = MPI.COMM_WORLD
		lb = dict()
		ub = dict()
		dual = dict()

		for i in range(1 , self.iNumAgents):
			lb[i] = comm.recv(source = i, tag = Constants.TREE_START_LB)
			ub[i] = comm.recv(source = i, tag = Constants.TREE_START_UB)
			dual[i] = comm.recv(source = i, tag = Constants.TREE_START_DUAL)				
		
		self.tree = SortedDict( [ (-float('inf') , Node( lb, ub , dual) ) ] )
			
	def get_next_node(self):
		return self.tree.popitem(last = False)

	def send_agent_node_info(self , node):
		comm = MPI.COMM_WORLD
		for i in range(1 , self.iNumAgents):
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
		for i in range(1 , self.iNumAgents):
			data = comm.recv(source = i , tag = Constants.GET_PRIMAL_SOLUTION)
			obj = obj + data[0];
			feasible = feasible * data[1]
			primal_sol[i] = data[2:]
			data = comm.recv(source = i , tag = Constants.GET_DUAL_VARS)
			dual_sol[i] = data
		return [obj , feasible, primal_sol, dual_sol]

	def create_new_nodes(self, obj, primal_sol , dual_sol , node):
		
		max = 0

		for i in range(1 , self.iNumAgents):	
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
		while (len(self.tree) > 0 ):
			[obj , node] = self.get_next_node()

			if(obj >= self.UB):
				continue

			self.send_agent_node_info(node)
			[obj , feasible, primal_sol, dual_sol] = self.get_solution_from_agents()
			self.create_new_nodes(obj, primal_sol , dual_sol , node)

			if (feasible) and (obj < self.UB):
				self.UB = obj
	

if __name__ == '__main__':
	
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	iNumAgents = 4;

	if (Constants.BB_TREE_ID == rank):
		tree = BB_Tree_Agent(iNumAgents)
		node = tree.get_next_node()
		print node[1].lb[1]
		print node[1].ub[2]
		#print node[1].dual[3]
		tree.begin_optimization()
	else:
		neigh = [rank]
		agent = Traffic_Agent(neigh, rank * 10)
		
		lb = []
		lb.append([])
		lb[0].append(rank)
		lb[0].append(rank)
		lb[0].append(rank)
		lb.append([])
		lb[1].append(2 * rank)
		lb[1].append(2 * rank)
		lb[1].append(2 * rank)

		ub = []
		ub.append([])
		ub[0].append(10*rank)
		ub[0].append(10*rank)
		ub[0].append(10*rank)
		ub.append([])
		ub[1].append(10*rank)
		ub[1].append(10*rank)
		ub[1].append(10*rank)

		dual = []
		dual.append([])
		dual[0].append(100 * rank)

		agent.set_new_data(lb, ub, dual)
		
		#agent.parse_file(%FileName)

		agent.construct_prob()
		agent.solve_problems()

		#comm.send([1, MPI.INT], dest = Constants.BB_TREE_ID , tag = Constants.READY_START)

	
	






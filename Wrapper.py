from sortedcontainers import SortedDict
from mpi4py import MPI
from cvxpy import *
import math
import Constants
import Solver
import numpy as np
import Network
from copy import deepcopy
import operator
import random

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
	
	def __init__(self, lb, ub , depth , pred_val):
		self.lb = deepcopy(lb)
		self.ub = deepcopy(ub)	
		self.depth = depth 
		self.pred_val = pred_val


def get_most_dev_index(dev , dev1, dev2, dev3):
		if(dev == dev1):
			return 1
		elif(dev == dev2):
			return 2
		else:
			return 3	

def get_ratio(num , primal_sol , lb , ub , mst_inf_ag , mst_inf_que):
	lb_val = lb[mst_inf_ag][mst_inf_que][num][0]
	ub_val = ub[mst_inf_ag][mst_inf_que][num][0]
	val = primal_sol[mst_inf_ag][mst_inf_que][num][0]

	if(ub_val == lb_val):
		return 0
	elif( ub_val - lb_val < 10 ** (-4) ):
		return 0	


	return ( (val - lb_val)/ (ub_val - lb_val)) * ( (ub_val - val)/(ub_val - lb_val) )

def get_vol_split(num1, num2 , mst_inf_ind ,primal_sol , lb , ub , mst_inf_ag , mst_inf_que):
	val1 = ub[mst_inf_ag][mst_inf_que][num1][0] - lb[mst_inf_ag][mst_inf_que][num1][0]
	val2 = ub[mst_inf_ag][mst_inf_que][num2][0] - lb[mst_inf_ag][mst_inf_que][num2][0]
	val3 = ub[mst_inf_ag][mst_inf_que][mst_inf_ind][0] - lb[mst_inf_ag][mst_inf_que][mst_inf_ind][0]

	val = max(val1, val2, val3)

	if(val == val1):
		return num1 , (ub[mst_inf_ag][mst_inf_que][num1][0] + lb[mst_inf_ag][mst_inf_que][num1][0])/2.0
	elif(val == val2):
		return num2	, (ub[mst_inf_ag][mst_inf_que][num2][0] + lb[mst_inf_ag][mst_inf_que][num2][0])/2.0
	else:
		return mst_inf_ind , (ub[mst_inf_ag][mst_inf_que][mst_inf_ind][0] + lb[mst_inf_ag][mst_inf_que][mst_inf_ind][0])/2.0

def get_best_split_ratio(num1 , num2 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind):
	ratio1 = get_ratio(num1 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que )
	ratio2 = get_ratio(num2 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que )
	ratio3 = get_ratio(mst_inf_ind , primal_sol , lb, ub , mst_inf_ag , mst_inf_que )

	ratio = max([ratio1 , ratio2, ratio3])

	if(ratio > 0.10):
		#print 'Non Volume Split'
		if(ratio == ratio1):
			return num1 , primal_sol[mst_inf_ag][mst_inf_que][num1][0]
		elif(ratio ==  ratio2):
			return num2 , primal_sol[mst_inf_ag][mst_inf_que][num2][0]
		else:
			return mst_inf_ind , primal_sol[mst_inf_ag][mst_inf_que][mst_inf_ind][0]		
	else:
		#print 'Volume Split'
		return get_vol_split(num1, num2 , mst_inf_ind ,primal_sol , lb , ub , mst_inf_ag , mst_inf_que)
				
def get_split_point(idx_1 , idx_2 , idx_3 , idx_4, idx_5, sp_indx , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind):
	if(idx_1):
		if(1 == sp_indx):
			return get_best_split_ratio(0 , 8 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)
		elif(2 == sp_indx):
			return get_best_split_ratio(13 , 5 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)
		else:
			return get_best_split_ratio(12 , 3 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)
	elif(idx_2):
		return get_best_split_ratio(5 , 10 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)	
	elif(idx_3):
		return get_best_split_ratio(3 , 7 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)	
	elif(idx_4):
		return get_best_split_ratio(9 , 3 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)	
	elif(idx_5):
		return get_best_split_ratio(11 , 4 , primal_sol , lb, ub , mst_inf_ag , mst_inf_que , mst_inf_ind)	


def set_data(dev , i , idx , j):
	return dev , i , idx, j

def set_truth(bInd1 , bInd2, bInd3, bInd4, bInd5):
	return 	bInd1 , bInd2, bInd3, bInd4, bInd5

class BB_Tree_Agent(object):
	
	def __init__(self , iNumAgents):
		self.iNumAgents = iNumAgents
		self.tree = [] 
		
	def get_node_info(self , agent_list):		
		lb = dict()
		ub = dict()
		
		for i in range(self.iNumAgents):
			lb[i] = agent_list[i].get_lb_bounds()
			ub[i] = agent_list[i].get_ub_bounds()			
		
		self.tree = SortedDict( [ (-float('inf') , Node( lb, ub , 0 , 0) ) ] )
	
	def get_next_node(self):
		return self.tree.popitem(last = False)

	def get_length(self):
		return len(self.tree)	

	def isEmpty(self):
		if(0 == len(self.tree)):
			return True
		else:	
			return False

	def create_new_nodes(self, opt, primal_sol , node , agent_list):		
		
		max_dev = 0
		ind_set = [ 1 , 6 , 0 , 20 , 21 ]

		best_split = 0
		idx_1 = False
		idx_2 = False
		idx_3 = False
		idx_4 = False
		idx_5 = False

		for i in range(self.iNumAgents):	
			for j in range(agent_list[i]._local_queue_num):
				for idx in ind_set:
					if(1 == idx):
						dev1 = np.absolute((primal_sol[i][j][0][0] * primal_sol[i][j][8][0]) - primal_sol[i][j][1][0])
						dev2 = np.absolute((primal_sol[i][j][13][0] * primal_sol[i][j][5][0]) - primal_sol[i][j][1][0])
						dev3 = np.absolute((primal_sol[i][j][12][0] * primal_sol[i][j][3][0]) - primal_sol[i][j][1][0])
						
						dev = max([dev1 , dev2 , dev3])	
						sp_indx = get_most_dev_index(dev , dev1, dev2, dev3)

						if(dev >= max_dev):
							max_dev , mst_inf_ag , mst_inf_ind , mst_inf_que = set_data(dev , i , idx, j)
							idx_1 , idx_2 , idx_3 , idx_4, idx_5 = set_truth(True, False, False, False, False)							

					elif(6 == idx):
						dev = np.absolute((primal_sol[i][j][5][0] * primal_sol[i][j][10][0]) - primal_sol[i][j][6][0])

						if(dev >= max_dev):
							max_dev , mst_inf_ag , mst_inf_ind , mst_inf_que = set_data(dev , i , idx, j)
							idx_1 , idx_2 , idx_3 , idx_4, idx_5 = set_truth(False, True, False, False, False)

					elif(0 == idx):
						dev = np.absolute((primal_sol[i][j][3][0] * primal_sol[i][j][7][0]) - primal_sol[i][j][0][0])	
						
						if(dev >= max_dev):
							max_dev , mst_inf_ag , mst_inf_ind , mst_inf_que = set_data(dev , i , idx, j)
							idx_1 , idx_2 , idx_3 , idx_4, idx_5 = set_truth(False, False, True, False, False)

					elif(20 == idx):	
						dev = np.absolute((primal_sol[i][j][9][0] * primal_sol[i][j][3][0]) - primal_sol[i][j][20][0])	
						
						if(dev >= max_dev):
							max_dev , mst_inf_ag , mst_inf_ind , mst_inf_que = set_data(dev , i , idx, j)
							idx_1 , idx_2 , idx_3 , idx_4, idx_5 = set_truth(False, False, False, True, False)	

					elif(21 == idx):

						dev = np.absolute((primal_sol[i][j][11][0] * primal_sol[i][j][4][0]) - primal_sol[i][j][21][0])	
						
						if(dev >= max_dev):
							max_dev , mst_inf_ag , mst_inf_ind , mst_inf_que = set_data(dev , i , idx, j)
							idx_1 , idx_2 , idx_3 , idx_4, idx_5 = set_truth(False, False, False, False, True)	


		if(max_dev <= 10^(-2)):
			return True

		split_index , split_val = get_split_point(idx_1 , idx_2 , idx_3 , idx_4, idx_5, sp_indx , primal_sol , node.lb, node.ub , mst_inf_ag , mst_inf_que , mst_inf_ind)

		print mst_inf_ag , mst_inf_que, split_index , node.lb[mst_inf_ag][mst_inf_que][split_index][0] , split_val, node.ub[mst_inf_ag][mst_inf_que][split_index][0]
		#print node.lb[mst_inf_ag][mst_inf_que][mst_inf_ind][0] , primal_sol[mst_inf_ag][mst_inf_que][mst_inf_ind][0] , node.ub[mst_inf_ag][mst_inf_que][mst_inf_ind][0]					

		if(split_val < node.lb[mst_inf_ag][mst_inf_que][split_index][0]):
			print 'Occurence1'
			split_val = node.lb[mst_inf_ag][mst_inf_que][split_index][0]
		elif(split_val > node.ub[mst_inf_ag][mst_inf_que][split_index][0]):
			print 'Occurence2'	
			split_val = node.ub[mst_inf_ag][mst_inf_que][split_index][0]

		ub1 = deepcopy(node.ub)
		ub1[mst_inf_ag][mst_inf_que][split_index][0] = split_val
			
		lb2 = deepcopy(node.lb)
		lb2[mst_inf_ag][mst_inf_que][split_index][0] = split_val

		key1 = 1000 * (opt + random.random() * Constants.OBJ_PERT_FACTOR)
		key2 = 1000 * (opt + random.random() * Constants.OBJ_PERT_FACTOR)

		self.tree[key1] = Node(lb2, node.ub , node.depth + 1 , opt)	
		self.tree[key2] = Node(node.lb, ub1 , node.depth + 1 , opt)

		#print key1 , node.lb[mst_inf_ag][mst_inf_que][split_index][0], ub1[mst_inf_ag][mst_inf_que][split_index][0]
		#print key2 , lb2[mst_inf_ag][mst_inf_que][split_index][0], node.ub[mst_inf_ag][mst_inf_que][split_index][0]
		
		return False
		
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
	
	#comm = MPI.COMM_WORLD
	#rank = comm.Get_rank()
	iNumAgents = 4;

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

	for i in range(0,4):
		t = agent_list[i]
		t.init_queue_solver_vars()
	
	tree = BB_Tree_Agent(4)
	tree.get_node_info(agent_list)

	best_sol = float('inf')

	while(tree.isEmpty() == False):

		constraints = []
		obj = 0	
		
		node_key = tree.get_next_node()
		
		key = node_key[0]
		node = node_key[1]

		if(key > best_sol):
			continue


		for i in range(0 , 4):		
			agent_list[i].set_lower_upper_bound(node)

		for i in range(0 , 4):
			t = agent_list[i]
			constraints += t.get_all_constraints()
			constraints += t.get_consensus_constraints()
			constraints += t.get_coupling_constraints(agent_list)
			obj = obj + t.get_primal_objective()
		
		opt = Minimize(obj)
		prob = Problem(opt , constraints)
		prob.solve(solver = ECOS , max_iters = 20000 , reltol = 10**(-9) , feastol = 10**(-6) , abstol = 10**(-7), verbose = False)
		
		#print 'Depth: ', node.depth , 'Opt Val: ', opt.value , 'Pred_val: ', node.pred_val, 'Accuracy: ' , prob.status , 'Tree Length: ', tree.get_length() 

		if('infeasible' == prob.status):
			print 'Depth: ', node.depth , 'Opt Val:', opt.value , 'Pred_val:', node.pred_val, 'Accuracy:' , prob.status , 'Tree Length:', tree.get_length()
			continue

		print 'Depth: ', node.depth , 'Delta Improv:', 1000 * (opt.value - node.pred_val), 'Accuracy:' , prob.status , 'Tree Length:', tree.get_length() 

		primal_sol = dict()
		for i in range(iNumAgents):
			primal_sol[i] = agent_list[i].get_sol()
		
		#Hsu-Chieh
		#call function, get the objective, if (objective < best_solution) best_solution = objective
		#hsu_chieh = ()

		# if(best_sol > hsu_chieh):
		# 	hsu_chieh = best_sol

		feas_sol = tree.create_new_nodes(opt.value , primal_sol , node , agent_list)	

		if(feas_sol):
			print 'Perfectly feasible solution'
			print 'Perfectly feasible solution'
			print 'Perfectly feasible solution'







		




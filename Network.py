import numpy as np
from mpi4py import MPI
from cvxpy import *
 

class TrafficNetwork(object):
	def __init__(self, agent_list, adjacent_matrix):
		self._vertex = agent_list
		self._vertex_num = len(agent_list)
		self._adjacent_matrix = adjacent_matrix

	def connect_queues(self):
		for vi in range(self._vertex_num):
			for vj in range(self._vertex_num):
				if self._adjacent_matrix[vi,vj] != ():
					self._vertex[vi].set_neighbors(self._vertex[vj], self._adjacent_matrix[vi,vj])					


class TrafficAgentModel(object):
	def __init__(self, agent_id, local_queue_num):
		self._agent_id = agent_id
		self._local_queue_num = local_queue_num
		self._local_queue = []
		self.init_local_queue()
		self._neighbor_agents = []

	def init_local_queue(self):
		for i in range(self._local_queue_num):
			queue = TrfficQueue(self._agent_id, i)
			self._local_queue.append(queue)

	def get_local_queue(self, queue_id):
		return self._local_queue[queue_id]


	def set_neighbors(self, neighbor, edge):
		
		if neighbor not in self._neighbor_agents:
			self._neighbor_agents.append(neighbor)

		for i in edge[0]:
			self._local_queue[edge[1]].set_upstream(neighbor.get_local_queue(edge[0][i]))
			neighbor.get_local_queue(i).set_downstream(self._local_queue[edge[1]])
	
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

	
class TrfficQueue(object):
	def __init__(self, agent_id, queue_id):
		self._agent_id = agent_id
		self._queue_id = queue_id

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


		'''constant'''
		self._up_turning_prop = []
		self._down_turning_prop = []
		self._ext_arr_rate = 0

	def init_queue(self):

	def set_upstream(self, queue):
		self._upstream_queue.append(queue)

	def set_downstream(self, queue):
		self._downstream_queue.append(queue)



def main():
	queue_num = 4
	t1 = TrafficAgentModel(1,queue_num)
	t2 = TrafficAgentModel(2,queue_num)
	t3 = TrafficAgentModel(3,queue_num)
	t4 = TrafficAgentModel(4,queue_num)
	agent_list = [t1,t2,t3,t4]
	adjacent_matrix = (((),(),(),())
	TrafficNetwork(agent_list, adjacent_matrix)



if __name__ == '__main__':
	main()




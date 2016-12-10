from cvxpy import *
import numpy as np
import Constants
from mpi4py import MPI

def Variable_Initialization(Up_queue, Down_queue):
	Vars = []
	Dual = dict()

	Vars.append([])
	Vars[0].append(Variable(1)) # 0 - \lambda_i
	
	Vars.append([])
	Vars[1].append(Variable(1)) # 1 - \lambda_i^{eff}
	
	Vars.append([])
	Vars[2].append(Variable(1)) # 2 - P(N_i = k_i)
	
	Vars.append([])
	Vars[3].append(Variable(1)) # 3 - \mu_i^{eff}
	
	Vars.append([])
	Vars[4].append(Variable(1)) # 4 - \mu_i

	Vars.append([])
	Vars[5].append(Variable(1)) # 5 - \tilde{\mu_i}

	Vars.append([])
	Vars[6].append(Variable(1)) # 6 - \mathcal{P}_i

	Vars.append([])
	Vars[7].append(Variable(1)) # 7 - \rho_i

	Vars.append([])
	Vars[8].append(Variable(1)) # 8 - v_i

	Vars.append([])
	Vars[9].append(Variable(1)) # 9 - d_i

	Vars.append([])
	Vars[10].append(Variable(1)) # 10 - t_i

	Vars.append([])
	Vars[11].append(Variable(1)) # 11 - w_i

	Vars.append([])
	Vars[12].append(Variable(1)) # 12 - f_i

	Vars.append([])
	Vars[13].append(Variable(1)) # 13 - r_i

	Vars.append(dict())
	Dual['A'] = dict()
	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		Vars[14][(agent_id, queue_id)] = 0			#14 - Az_i^j 
		Dual['A'][(agent_id, queue_id)] = 0 

	Vars.append([])	
	Vars[15].append(0)								#15 - Bz_i
	Dual['B'] = 0 

	Vars.append([])	
	Vars[16].append(0)								#16 - Cz_i
	Dual['C'] = 0

	Vars.append(dict())
	Dual['D'] = dict()
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		Vars[17][(agent_id, queue_id)] = 0			#17 - Dz_i^j
		Dual['D'][(agent_id, queue_id)] = 0

	Vars.append([])			
	Vars[18].append(0)								#18 - EZ_i
	Dual['E'] = 0
	
	Vars.append(dict())
	Dual['F'] = dict()
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		Vars[19][(agent_id, queue_id)] = 0			#19 - Fz_i^j 
		Dual['F'][(agent_id, queue_id)] = 0
	 
	Dual['1'] = 0
	Dual['4'] = 0
	Dual['5'] = 0	

	return Vars , Dual

def Create_constraints_for_queue(Vars, Up_queue, Down_queue , lb, ub , queue_cap):
	constraints = []
	
	eqn_bounds = create_lower_upper_bound_constraints(Vars , lb , ub)
	constraints += eqn_bounds

	eqn1 = Create_eqn_1(Vars, lb, ub)
	constraints += eqn1

	eqn3 = Create_eqn_3(Vars, lb, ub)
	constraints += eqn3

	eqn4 = Create_eqn_4_partial(Vars, lb, ub)
	constraints += eqn4

	eqn6 = Create_eqn_6(Vars, lb, ub , queue_cap)
	constraints += eqn6

	eqn7 = Create_eqn_7(Vars, lb, ub)
	constraints += eqn7

	return constraints

def create_lower_upper_bound_constraints(Vars , lb , ub):
	eqns = []
	for i in range(0 , 14):
		eqns.append( Vars[i][0] >= lb[i][0] )
		eqns.append( Vars[i][0] <= ub[i][0] )
	return eqns	

def Create_eqn_1(Vars , lb, ub):
	return Construct_McCormick( Vars[0][0] , Vars[8][0] , Vars[1][0], lb[0][0], ub[0][0] , lb[8][0] , ub[8][0] )

def Create_eqn_3(Vars , lb, ub):
	eqns = []

	eqns.append( Vars[11][0] + Vars[10][0] ==  Vars[9][0] )     #eqn above 9
	eqns += Construct_Fractional_Envelope( Vars[9][0] , Vars[3][0] , lb[9][0] , ub[9][0] , lb[3][0] , ub[3][0] )        #eqn 9
	eqns += Construct_Fractional_Envelope( Vars[11][0] , Vars[4][0] , lb[11][0] , ub[11][0] , lb[4][0] , ub[4][0] ) 	  #eqn10	
	eqns += Construct_McCormick( Vars[5][0] , Vars[10][0] , Vars[6][0] , lb[5][0] , ub[5][0] , lb[10][0] , ub[10][0] )   #eqn 11

	return eqns

def Create_eqn_4_partial(Vars , lb, ub):
	eqns = []
	eqns += Construct_McCormick( Vars[13][0] , Vars[5][0] , Vars[1][0] , lb[13][0] , ub[13][0] , lb[5][0] , ub[5][0] )   #eqn 12
	eqns += Construct_McCormick( Vars[12][0] , Vars[3][0] , Vars[1][0] , lb[12][0] , ub[12][0] , lb[3][0] , ub[3][0] )   #eqn 13
	return eqns

def Create_eqn_6(Vars , lb , ub , queue_cap ):
	eqns = []

	var_lb = lb[7][0]
	var_ub = ub[7][0]

	K1 =  ((1 - var_lb)/var_lb) * ( (1/(1 - (var_lb ** (queue_cap+1)))) - 1 )  
	K2 =  ((1 - var_ub)/var_ub) * ( (1/(1 - (var_ub ** (queue_cap+1)))) - 1 )

	eqns.append( ((K2 - K1)*( Vars[7][0] - var_lb )) + ((var_ub - var_lb)*K1) >=  (var_ub - var_lb)*Vars[2][0] )

	#can be changed to improve accuracy
	eqns.append( Vars[2][0] >= K1 )

	return eqns

def Create_eqn_7(Vars, lb, ub):
	return Construct_McCormick(Vars[7][0] , Vars[3][0] , Vars[0][0] , lb[7][0] , ub[7][0] , lb[3][0] , ub[3][0] )

#  var1 * var2 = var3
def Construct_McCormick(var1 , var2 , var3, lb1 , ub1 , lb2 , ub2 ):
	eqns = []
	eqns.append( ( lb1 * var2 ) + ( lb2 * var1 ) - ( lb1 * lb2 ) <= var3 )
	eqns.append( ( ub1 * var2 ) + ( ub2 * var1 ) - ( ub1 * ub2 ) <= var3 )
	eqns.append( ( ub1 * var2 ) + ( lb2 * var1 ) - ( ub1 * lb2 ) >= var3 )
	eqns.append( ( lb1 * var2 ) + ( ub2 * var1 ) - ( lb1 * ub2 ) >= var3 )

	return eqns	   	

#  var1 * var2 = 1
def Construct_Fractional_Envelope(var1 , var2 , lb1 , ub1 , lb2 , ub2 ):
	eqns = []
	eqns.append( (var1 * (ub2 - lb2)) + (var2 * (ub1 - lb1)) <= (ub2 * ub1) - (lb1 * ub1) )
	eqns.append( (ub2 * lb1 * lb1 ) + (lb1 - var1) <= (var2 * lb1 * lb1) )
	eqns.append( (lb2 * ub1 * ub1 ) + (ub1 - var1) <= (var2 * ub1 * ub1) )

	return eqns

def Get_total_queue_objective(Vars, Duals, Up_queue, Down_queue , ext_arr_rate , turn_prop , up_turn_prop):
	obj = inv_pos(1 - Vars[7][0]) - 1
	obj = obj + Get_dual_objective_eqn_2(Vars, Duals, Down_queue , ext_arr_rate , turn_prop)
	obj = obj + Get_dual_objective_eqn_4(Vars, Duals, Up_queue)
	obj = obj + Get_dual_objective_eqn_5(Vars, Duals, up_turn_prop , Up_queue)
	return obj

def Get_dual_objective_eqn_2(Vars, Duals , Down_queue , ext_arr_rate , turn_prop):	
	obj = Get_dual_objective_eqn_A(Vars , Duals , Down_queue , turn_prop)
	obj = obj + Get_dual_objective_eqn_B(Vars , Duals , ext_arr_rate)
	return obj	

def Get_dual_objective_eqn_4(Vars , Duals , Up_queue):
	obj  = 	Get_dual_objective_eqn_C(Vars , Duals)
	obj = obj + Get_dual_objective_eqn_D(Vars, Duals , Up_queue)
	return obj	

def Get_dual_objective_eqn_5(Vars, Duals, turn_prop , Up_queue):
	obj = Get_dual_objective_eqn_E(Vars, Duals)
	obj = obj + Get_dual_objective_eqn_F(Vars, Duals, turn_prop , Up_queue)
	return obj		

def Get_dual_objective_eqn_A(Vars , Duals , Down_queue , turn_prop):	
	obj = 0
	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		var = (turn_prop[(agent_id, queue_id)] * Vars[1][0] ) - Vars[14][(agent_id, queue_id)]
		obj = obj + Duals['A'][(agent_id, queue_id)] * var 
		obj = obj + Constants.ADMM_PEN * square( var )/2.0
	return obj	

def Get_dual_objective_eqn_B(Vars , Duals , arr_rate):	
	var = Vars[1][0] - (arr_rate * Vars[8][0]) - Vars[15][0]	
	obj = ( Duals['B'] * var ) + Constants.ADMM_PEN * square(var)/2.0 
	return obj

def Get_dual_objective_eqn_C(Vars, Duals):
	var = Vars[13][0] - Vars[16][0]
	obj = ( Duals['C'] * var ) + Constants.ADMM_PEN * square(var)/2.0
	return obj	

def Get_dual_objective_eqn_D(Vars, Duals , Up_queue):
	obj = 0
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Vars[12][0] - Vars[17][(agent_id , queue_id)]
		obj = obj + Duals['D'][(agent_id, queue_id)] * var
		obj = obj + Constants.ADMM_PEN * square(var)/2.0 
	return obj		

def Get_dual_objective_eqn_E(Vars, Duals):	
	var = Vars[6][0] - Vars[18][0]
	obj = (Duals['E'] * var) + Constants.ADMM_PEN * square(var)/2.0 
	return obj

def Get_dual_objective_eqn_F(Vars, Duals, turn_prop , Up_queue):
	obj = 0
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Vars[19][(agent_id , queue_id)] - (turn_prop[(agent_id , queue_id)] * Vars[2][0])	
		obj = obj + Duals['F'][(agent_id , queue_id)] * var
		obj = obj + Constants.ADMM_PEN * square(var)/2.0 
	return obj

def Update_Dual_Vars(Vars , Duals, Up_queue, Down_queue , ext_arr_rate , turn_prop , up_turn_prop):
	Update_Dual_Vars_eqn_2(Vars , Duals, Down_queue , ext_arr_rate , turn_prop)	
	Update_Dual_Vars_eqn_4(Vars, Duals, Up_queue)
	Update_Dual_Vars_eqn_5(Vars, Duals, up_turn_prop)

def Update_Dual_Vars_eqn_2(Vars , Duals, Down_queue , ext_arr_rate , turn_prop):
	Update_Dual_Vars_eqn_A(Vars , Duals , Down_queue , turn_prop)	
	Update_Dual_Vars_eqn_B(Vars , Duals , ext_arr_rate)

def Update_Dual_Vars_eqn_4(Vars, Duals, Up_queue):	
	Update_Dual_Vars_eqn_C(Vars , Duals)
	Update_Dual_Vars_eqn_D(Vars, Duals , Up_queue)

def Update_Dual_Vars_eqn_5(Vars, Duals, turn_prop):
	Update_Dual_Vars_eqn_E(Vars, Duals)
	Update_Dual_Vars_eqn_F(Vars, Duals, turn_prop)

def Update_Dual_Vars_eqn_A(Vars , Duals , Down_queue , turn_prop):
	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		var = (turn_prop[(agent_id, queue_id)] * Vars[1][0] ) - Vars[14][(agent_id, queue_id)]
		Dual['A'][(agent_id, queue_id)]  = Dual['A'][(agent_id, queue_id)] + ( Constants.ADMM_PEN * var )
		
def Update_Dual_Vars_eqn_B(Vars , Duals , arr_rate):	
	var = Vars[1][0] - (arr_rate * Vars[8][0]) - Vars[15][0]	
	Dual['B'] = Dual['B'] + ( Constants.ADMM_PEN * var ) 
		
def Update_Dual_Vars_eqn_C(Vars, Duals):
	var = Vars[13][0] - Vars[16][0]
	Duals['C'] = Duals['C'] + ( Constants.ADMM_PEN * var ) 

def Update_Dual_Vars_eqn_D(Vars, Duals , Up_queue):
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Vars[12][0] - Vars[17][(agent_id , queue_id)]
		Dual['D'][(agent_id, queue_id)] = Dual['D'][(agent_id, queue_id)] + ( Constants.ADMM_PEN * var )

def Update_Dual_Vars_eqn_E(Vars, Duals):	
	var = Vars[6][0] - Vars[18][0]
	Duals['E'] = Duals['E'] + (Constants.ADMM_PEN * var) 

def Update_Dual_Vars_eqn_F(Vars, Duals, turn_prop):
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Vars[19][(agent_id , queue_id)] - (turn_prop[(agent_id , queue_id)] * Vars[2][0])	
		Dual['F'][(agent_id , queue_id)] = Dual['F'][(agent_id , queue_id)] + ( Constants.ADMM_PEN * var )	
	
def send_rel_vars(Vars , Duals, Up_queue , Down_queue , turn_prop , My_queue_id , lb, ub , up_turn_prop):
	send_rel_vars_eqn_A( Vars , Duals , Down_queue , turn_prop , My_queue_id , lb , ub )
	send_rel_vars_eqn_D(Vars , Duals , Up_queue , My_queue_id ,lb , ub)
	send_rel_vars_eqn_F(Vars , Duals , Up_queue , up_turn_prop , My_queue_id , lb , ub )

def send_rel_vars_eqn_A( Vars , Duals , Down_queue , turn_prop , My_queue_id , lb , ub):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		TAG = str(rank)+'_'+str(queue_id)+'_' +str(My_queue_id)+'A'
		data = [ turn_prop[(agent_id , queue_id)] * Vars[1][0].value , Duals['A'][(agent_id, queue_id)] , lb[14][(agent_id , queue_id)], ub[14][(agent_id , queue_id)] ]
		comm.send( data , dest = agent_id , tag = TAG )	

def send_rel_vars_eqn_D(Vars , Duals , Up_queue , My_queue_id ):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()		
	
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		TAG = str(rank) + '_' + str(queue_id) +'_' +str(My_queue_id)+'D'
		data = [ Vars[12][0].value , Dual['D'][(agent_id, queue_id)] , lb[17][(agent_id, queue_id)], ub[17][(agent_id, queue_id)] ]
		comm.send( data , dest = agent_id , tag = TAG )

def send_rel_vars_eqn_F(Vars , Duals , Up_queue , turn_prop , My_queue_id ):		
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		TAG = str(rank) + '_' + str(queue_id) + '_' +str(My_queue_id) + 'F'
		data = [ turn_prop[(agent_id , queue_id)] * Vars[2][0].value , Dual['F'][(agent_id, queue_id)], lb[19][(agent_id, queue_id)], ub[19][(agent_id, queue_id)] ]
		comm.send( data , dest = agent_id , tag = TAG )	

def receive_rel_vars(Vars , Duals, Up_queue , Down_queue , My_queue_id ):
	neigh_vars_data = dict()
	receive_rel_vars_A(Vars , Duals , Up_queue , My_queue_id , neigh_vars_data)
	receive_rel_vars_D(Vars , Duals , Down_queue , My_queue_id , neigh_vars_data)
	receive_rel_vars_F(Vars , Duals , Down_queue , My_queue_id , neigh_vars_data)
	return neigh_vars_data

def receive_rel_vars_A(Vars , Duals , Up_queue , My_queue_id , neigh_vars_data):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	neigh_vars_data['A'] = dict()

	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		TAG = str(agent_id)+'_'+str(My_queue_id)+ '_' + str(queue_id) + 'A'
		data = comm.recv(source = agent_id , tag = TAG)	
		neigh_vars_data['A'][(agent_id , queue_id)] = data

def receive_rel_vars_D(Vars , Duals , Down_queue , My_queue_id , neigh_vars_data):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	neigh_vars_data['D'] = dict()

	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		TAG = str(agent_id)+'_'+str(My_queue_id)+ '_' + str(queue_id) + 'D'
		data = comm.recv(source = agent_id , tag = TAG)	
		neigh_vars_data['D'][(agent_id , queue_id)] = data

def receive_rel_vars_F(Vars , Duals , Down_queue , My_queue_id , neigh_vars_data):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	neigh_vars_data['F'] = dict()

	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		TAG = str(agent_id)+'_'+str(My_queue_id)+ '_' + str(queue_id) + 'F'
		data = comm.recv(source = agent_id , tag = TAG)	
		neigh_vars_data['F'][(agent_id , queue_id)] = data		

def solve_coupling_eqns_send_sols(Vars , Duals , Up_queue , Down_queue , My_queue_id , neigh_data , arr_rate ):
	solve_send_eqn_2( Vars , Duals , Up_queue , My_queue_id , neigh_data , arr_rate )
	solve_send_eqn_4( Vars , Duals , Down_queue , My_queue_id , neigh_data )
	solve_send_eqn_5( Vars , Duals , Down_queue , My_queue_id , neigh_data )

def solve_send_eqn_2( Vars , Duals , Up_queue , My_queue_id , neigh_data , arr_rate ):
	Conses_neigh_vars = dict()
	solve_eqn_2	( Vars , Duals , Up_queue , My_queue_id , neigh_data , arr_rate , Conses_neigh_vars )
	send_solved_eqn(Conses_neigh_vars , Up_queue , My_queue_id , 2)	

def solve_send_eqn_4( Vars , Duals , Down_queue , My_queue_id , neigh_data ):
	Conses_neigh_vars = dict()
	solve_eqn_4( Vars , Duals , Down_queue , My_queue_id , neigh_data , Conses_neigh_vars )
	send_solved_eqn(Conses_neigh_vars , Down_queue , My_queue_id , 4)

def solve_send_eqn_5( Vars , Duals , Down_queue , My_queue_id , neigh_data ):
	Conses_neigh_vars = dict()
	solve_eqn_5( Vars , Duals , Down_queue , My_queue_id , neigh_data , Conses_neigh_vars )
	send_solved_eqn(Conses_neigh_vars , Down_queue , My_queue_id , 5)	

def solve_eqn_2( Vars , Duals , Up_queue , My_queue_id , neigh_data , arr_rate , Conses_neigh_vars ):
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		Conses_neigh_vars[(agent_id , queue_id)] = Variable(1)

	constraints = []

	BZ_Var = Variable(1)
	var = Vars[1][0] - (arr_rate * Vars[8][0]) - BZ_Var	
	obj = ( Duals['B'] * var ) + ( Constants.ADMM_PEN * var/2.0 )

	constraints.append(BZ_Var >= lb[15][0])
	constraints.append(BZ_Var <= ub[15][0])

	coup_var = -1.0 * BZ_Var

	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Conses_neigh_vars[(agent_id , queue_id)]
		agent_obj  = neigh_data['A'][(agent_id , queue_id)][1] * (neigh_data['A'][(agent_id , queue_id)][0] - var)
		agent_obj = agent_obj + Constants.ADMM_PEN * square(neigh_data['A'][(agent_id , queue_id)][0] - var)/2.0

		constraints.append( var >= neigh_data['A'][(agent_id , queue_id)][2] )
		constraints.append( var <= neigh_data['A'][(agent_id , queue_id)][3] )

		obj = obj + agent_obj
		coup_var = coup_var + var

	obj = obj + ( Duals['1'] * coup_var ) + Constants.ADMM_PEN * square(coup_var)/2.0
	prob = Minimize(obj , constraints)

	Vars[15][0] = BZ_Var.value
	#dual update of coupling constraint 1
	Duals['1'] = Duals['1'] + Constants.ADMM_PEN * ( coup_var.value )

def solve_eqn_4( Vars , Duals , Down_queue , My_queue_id , neigh_data , Conses_neigh_vars ):
	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		Conses_neigh_vars[(agent_id , queue_id)] = Variable(1)

	constraints = []
		
	CZ_Var = Variable(1)
	var = Vars[13][0] - CZ_Var
	obj = (Duals['C'] * var) + Constants.ADMM_PEN * square(var)/2.0

	constraints.append(BZ_Var >= lb[16][0])
	constraints.append(BZ_Var <= ub[16][0])

	coup_var = -1.0 * CZ_Var

	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		var = Conses_neigh_vars[(agent_id , queue_id)]
		agent_obj = neigh_data['D'][(agent_id , queue_id)][1] * (neigh_data['D'][(agent_id , queue_id)][0] - var)
		agent_obj = agent_obj + Constants.ADMM_PEN * square(neigh_data['D'][(agent_id , queue_id)][0] - var)/2.0

		constraints.append( var >= neigh_data['D'][(agent_id , queue_id)][2] )
		constraints.append( var <= neigh_data['D'][(agent_id , queue_id)][3] )

		obj = obj + agent_obj
		coup_var = coup_var + var

	obj = obj + (Duals['4'] * coup_var) + Constants.ADMM_PEN * square(coup_var)/2.0
	prob = Minimize(obj , constraints)
	Vars[16][0] = CZ_Var.value
	#dual update of coupling constraint 4
	Duals['4'] = Duals['4'] + Constants.ADMM_PEN * ( coup_var.value )	

def solve_eqn_5( Vars , Duals , Down_queue , My_queue_id , neigh_data , Conses_neigh_vars ):
	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		Conses_neigh_vars[(agent_id , queue_id)] = Variable(1)

	constraints = []	

	EZ_Var = Variable(1)
	var = Vars[6][0] - EZ_Var
	obj = (Duals['E'] * var) + Constants.ADMM_PEN * square(var)/2.0

	constraints.append(EZ_Var >= lb[18][0])
	constraints.append(EZ_Var <= ub[18][0])

	coup_var = -1.0 * EZ_Var	

	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		var = Conses_neigh_vars[(agent_id , queue_id)]
		agent_obj = neigh_data['F'][(agent_id , queue_id)][1] * ( var - neigh_data['F'][(agent_id , queue_id)][0] )
		agent_obj = agent_obj + Constants.ADMM_PEN * square( var - neigh_data['F'][(agent_id , queue_id)][0] )/2.0

		constraints.append(var >= neigh_data['F'][(agent_id , queue_id)][2])
		constraints.append(var <= neigh_data['F'][(agent_id , queue_id)][3])

		obj = obj + agent_obj
		coup_var = coup_var + var

	obj = obj + (Duals['5'] * coup_var) + Constants.ADMM_PEN * square(coup_var)/2.0
	prob = Minimize(obj , constraints)
	Vars[18][0] = EZ_Var.value	
	#dual update of coupling constraint 5
	Duals['5'] = Duals['5'] + Constants.ADMM_PEN * ( coup_var.value )	

def send_solved_eqn(Conses_neigh_vars , Queue , My_queue_id , eqn_type):
	comm = MPI.COMM_WORLD
	
	for i in range(len(Queue)):
		agent_id = Queue[i]._agent_id
		queue_id = Queue[i]._queue_id
		TAG = str(agent_id) + '_' + str(queue_id) + '_' + str(My_queue_id) + str(eqn_type)
		comm.send( Conses_neigh_vars[(agent_id , queue_id)].value , dest = agent_id , tag = TAG )

def recv_update_solved_vars(Vars , My_agent_id , My_queue_id , Up_queue , Down_queue):
	recv_updated_solved_vars(Vars , My_agent_id , My_queue_id , 2 , Down_queue , 14)
	recv_updated_solved_vars(Vars , My_agent_id , My_queue_id , 4 , Up_queue , 17)
	recv_updated_solved_vars(Vars , My_agent_id , My_queue_id , 5 , Up_queue , 19)

def recv_updated_solved_vars(Vars , My_agent_id , My_queue_id , eqn_type , Queue , var_index):
	comm = MPI.COMM_WORLD

	for i in range(len(Queue)):
		agent_id = Queue[i]._agent_id
		queue_id = Queue[i]._queue_id
		TAG = str(My_agent_id)+'_'+str(My_queue_id)+'_'+str(queue_id)+str(eqn_type)
		data = comm.recv(source = agent_id , tag = TAG)
		Vars[var_index][(agent_id, queue_id)] = data
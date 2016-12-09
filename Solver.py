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
	  
	Vars.append([])
	Vars[20].append(0)

	Vars.append([])
	Vars[21].append(0)

	return Vars , Dual

def Create_constraints_for_queue(Vars, Up_queue, Down_queue , lb, ub):
	constraints = []
	
	eqn1 = Create_eqn_1(Vars, lb, ub)
	constraints.append(eqn1)

	eqn3 = Create_eqn_3(Vars, lb, ub)
	constraints.append(eqn3)

	eqn4 = Create_eqn_4_partial(Vars, lb, ub)
	constraints.append(eqn4)

	eqn6 = Create_eqn_6(Vars, lb, ub)
	constraints.append(eqn6)

	eqn7 = Create_eqn_7(Vars, lb, ub)
	constraints.append(eqn7)

def Create_eqn_1(Vars , lb, ub):
	return Construct_McCormick( Vars[0][0] , Vars[0][8] , Vars[1][0], lb[0][0], ub[0][0] , lb[0][8] , ub[0][8] )

def Create_eqn_3(Vars , lb, ub):
	eqns = []

	eqns.append( Vars[11][0] + Vars[10][0] ==  Vars[9][0] )     #eqn above 9
	eqns.append( Construct_Fractional_Envelope( Vars[9][0] , Vars[3][0] , lb[9][0] , ub[9][0] , lb[3][0] , ub[3][0] ) )       #eqn 9
	eqns.append( Construct_Fractional_Envelope( Vars[11][0] , Vars[4][0] , lb[11][0] , ub[11][0] , lb[4][0] , ub[4][0] ) )	  #eqn10	
	eqns.append( Construct_McCormick( Vars[5][0] , Vars[10][0] , Vars[6][0] , lb[5][0] , ub[5][0] , lb[10][0] , ub[10][0] ) )  #eqn 11

	return eqns

def Create_eqn_4_partial(Vars , lb, ub):
	eqns = []

	eqns.append( Construct_McCormick( Vars[13][0] , Vars[5][0] , Vars[1][0] , lb[13][0] , ub[13][0] , lb[5][0] , ub[5][0] ) )  #eqn 12
	eqns.append( Construct_McCormick( Vars[12][0] , Vars[3][0] , Vars[1][0] , lb[12][0] , ub[12][0] , lb[3][0] , ub[3][0] ) )  #eqn 13
	return eqns

def Create_eqn_6(Vars , lb , ub ):
	eqns = []

	var_lb = lb[7][0]
	var_ub = ub[7][0]

	K1 =  ((1 - var_lb)/var_lb) * ( (1/(1 - (var_lb ** (Constants.QUEUE_LENGTH+1)))) - 1 )  
	K2 =  ((1 - var_ub)/var_ub) * ( (1/(1 - (var_ub ** (Constants.QUEUE_LENGTH+1)))) - 1 )

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

def Get_queue_objective(Vars, Duals, Up_queue, Down_queue , ext_arr_rate , turn_prop):
	obj = inv_pos(1 - Vars[7][0]) - 1
	obj = obj + Get_dual_objective_eqn_2(Vars, Duals, Down_queue , ext_arr_rate , turn_prop)
	obj = obj + Get_dual_objective_eqn_4(Vars, Duals, Up_queue)
	obj = obj + Get_dual_objective_eqn_5(Vars, Duals, turn_prop)
	return obj

def Get_dual_objective_eqn_2(Vars, Duals , Down_queue , ext_arr_rate , turn_prop):	
	obj = Get_dual_objective_eqnA(Vars , Duals , Down_queue , turn_prop)
	obj = obj + Get_dual_objective_eqnB(Vars , Duals , ext_arr_rate)
	return obj	

def Get_dual_objective_eqn_4(Vars , Duals , Up_queue):
	obj  = 	Get_dual_objective_eqn_C(Vars , Duals)
	obj = obj + Get_dual_objective_eqn_D(Vars, Duals , Up_queue)
	return obj	

def Get_dual_objective_eqn_5(Vars, Duals, turn_prop):
	obj = Get_dual_objective_eqn_E(Vars, Duals)
	obj = obj + Get_dual_objective_eqn_F(Vars, Duals, turn_prop)
	return obj		

def Get_dual_objective_eqn_A(Vars , Duals , Down_queue , turn_prop):	
	obj = norm(0)
	for i in range(len(Down_queue)):
		agent_id = Down_queue[i]._agent_id
		queue_id = Down_queue[i]._queue_id
		var = turn_prop[(agent_id, queue_id)] * Vars[1][0] ) - Vars[14][(agent_id, queue_id)]
		obj = obj + Dual['A'][(agent_id, queue_id)] * var 
		obj = obj + (Constants.ADMM_PEN * square( var )/2.0
	return obj	

def Get_dual_objective_eqn_B(Vars , Duals , arr_rate):	
	var = Vars[1][0] - (arr_rate * Vars[8][0]) - Vars[15][0]	
	obj = ( Dual['B'] * var ) + Constants.ADMM_PEN * square(var)/2.0 
	return obj

def Get_dual_objective_eqn_C(Vars, Duals):
	var = Vars[13][0] - Vars[16][0]
	obj = ( Duals['C'] * var ) + Constants.ADMM_PEN * square(var)/2.0
	return obj	

def Get_dual_objective_eqn_D(Vars, Duals , Up_queue):
	obj = norm(0)
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Vars[12][0] - Vars[17][(agent_id , queue_id)]
		obj = obj + Dual['D'][(agent_id, queue_id)] * var
		obj = obj + Constants.ADMM_PEN * square(var)/2.0 
	return obj		

def Get_dual_objective_eqn_E(Vars, Duals):	
	var = Vars[6][0] - Vars[18][0]
	obj = (Duals['E'] * var) + Constants.ADMM_PEN * square(var)/2.0 
	return obj

def Get_dual_objective_eqn_F(Vars, Duals, turn_prop):
	obj = norm(0)
	for i in range(len(Up_queue)):
		agent_id = Up_queue[i]._agent_id
		queue_id = Up_queue[i]._queue_id
		var = Vars[19][(agent_id , queue_id)] - (turn_prop((agent_id , queue_id)) * Vars[2][0])	
		obj = obj + Dual['F'][(agent_id , queue_id)] * var
		obj = obj + Constants.ADMM_PEN * square(var)/2.0 
	return obj

def Update_Dual_Vars(Vars , Duals, Up_queue, Down_queue , ext_arr_rate , turn_prop):
	Update_Dual_Vars_eqn_2(Vars , Duals, Down_queue , ext_arr_rate , turn_prop)	
	Update_Dual_Vars_eqn_4(Vars, Duals, Up_queue)
	Update_Dual_Vars_eqn_5(Vars, Duals, turn_prop)

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
		var = turn_prop[(agent_id, queue_id)] * Vars[1][0] ) - Vars[14][(agent_id, queue_id)]
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
		var = Vars[19][(agent_id , queue_id)] - (turn_prop((agent_id , queue_id)) * Vars[2][0])	
		Dual['F'][(agent_id , queue_id)] = Dual['F'][(agent_id , queue_id)] + ( Constants.ADMM_PEN * var )	
	
def send_rel_vars(Vars , Duals, Up_queue , Down_queue , turn_prop):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

def send_rel_vars_eqn_A( Vars , Duals ,  )





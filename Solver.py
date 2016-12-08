from cvxpy import *
import numpy as np
import Constants

def Variable_Initialization(self):
	Vars = []

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

	Vars.append([])
	for i in range(0 , len(self.upstream)):
		Vars[14].append(Variable(1))				#14 - Az_i^j 

	Vars.append([])	
	Vars[15].append(Variable(1))					#15 - Bz_i

	Vars.append([])	
	Vars[16].append(Variable(1))					#16 - Cz_i

	Vars.append([])
	for i in range(0 , len(self.downstream)):
		Vars[17].append(Variable(1))					#17 - Dz_i^j

	Vars.append([])			
	Vars[18].append(Variable(1))					#18 - EZ_i
	
	Vars.append([])
	for i in range(0 , len(self.upstream)):
		Vars[19].append(Variable(1))				#19 - Fz_i^j 
	  
	return Vars


def Create_constraints_for_queue(Vars, Up_queue, down_queue , lb, ub):
	constraints = []
	
	eqn1 = Create_eqn_1(Vars, lb, ub)
	constraints.append(eqn1)

	eqn3 = Create_eqn_3(Vars, lb, ub)
	constraints.append(eqn3)

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


from cvxpy import *
import numpy as np

def Variable_Initialization(self, iNumVars):
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

def Create_Constraint_Matrix(lb , ub, iNumEqns , iNumVars , Vars)
	
	constraints = [ (lb[0][0]*Vars[8][0]) + (lb[8][0]*Vars[0][0]) - (lb[0][0]*lb[8][0]) <= Vars[2][0],
					(ub[0][0]*Vars[8][0]) + (ub[8][0]*Vars[0][0]) - (ub[0][0]*ub[0][0]) <= Vars[2][0],	
					(ub[0][0]*Vars[8][0]) + (lb[8][0]*Vars[0][0]) - (ub[0][0]*lb[8][0]) >= Vars[2][0],
					(lb[0][0]*Vars[8][0]) + (ub[8][0]*Vars[0][0]) - (lb[0][0]*ub[8][0]) >= Vars[2][0],
				  
					
				  ] 


import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def damagePlastic1D(params, eps, dt):
	"""
	Input: Parameter array, strain array
	Output: stress array
	"""

	# Parameters
	E 			= params[0]
	C 			= params[1]
	gamma 		= params[2]
	s 			= params[3]
	sig0 		= params[4]
	tau0 		= params[5]
	epsDot0 	= params[6]
	m 			= params[7]

	# Exponential damage law 
	def d(epsP):
		return 1-np.exp(-s*abs(epsP))

	# Define stress
	def sigma(eps_i,epsP_i):
		return (1-d(epsP_i))*E*(eps_i-epsP_i)

	# Initialize
	nSteps = len(eps)
	sig = np.zeros(nSteps)
	epsPArray = np.zeros(nSteps)
	dmg = np.zeros(nSteps)
	epsPOld = 0

	# Integrate
	for i in range(nSteps):
		 
		# Define effective potential
		def f(epsP):
			term1 = 0.5*(1-d(epsP))*E*(eps[i]-epsP)*(eps[i]-epsP)
			term2 = C/(gamma+1)*pow(np.abs(epsP),gamma+1)
			term3 = sig0 * abs(epsP-epsPOld)
			term4 = tau0/(m+1) * pow(abs(epsP-epsPOld),m+1)/pow(dt*epsDot0,m)
			return term1 + term2 + term3 + term4

		# Solve for new plastic strain
		res = minimize(f,x0=epsPOld, method='Nelder-Mead',
			tol=1e-8, options={'initial_simplex':[[0],[1]]})
		epsP = res['x']
		sig[i] = sigma(eps[i],epsP)
		dmg[i] = d(epsP)
		epsPArray[i] = epsP
		epsPOld = epsP

	return sig,dmg,epsPArray

def damagePlastic1D(params, eps, dt):
	"""
	Input: Parameter array, strain array
	Output: stress array
	"""

	# Parameters
	E 			= params[0]
	C 			= params[1]
	gamma 		= params[2]
	s 			= params[3]
	sig0 		= params[4]
	tau0 		= params[5]
	epsDot0 	= params[6]
	m 			= params[7]

	# Exponential damage law 
	def d(epsP):
		return 1-np.exp(-s*abs(epsP))

	# Define stress
	def sigma(eps_i,epsP_i):
		return (1-d(epsP_i))*E*(eps_i-epsP_i)

	# Initialize
	nSteps = len(eps)
	sig = np.zeros(nSteps)
	epsPArray = np.zeros(nSteps)
	dmg = np.zeros(nSteps)
	epsPOld = 0

	# Integrate
	for i in range(nSteps):
		 
		# Define effective potential
		def f(epsP):
			term1 = 0.5*(1-d(epsP))*E*(eps[i]-epsP)*(eps[i]-epsP)
			term2 = C/(gamma+1)*pow(np.abs(epsP),gamma+1)
			term3 = sig0 * abs(epsP-epsPOld)
			term4 = tau0/(m+1) * pow(abs(epsP-epsPOld),m+1)/pow(dt[i]*epsDot0,m)
			return term1 + term2 + term3 + term4

		# Solve for new plastic strain
		res = minimize(f,x0=epsPOld, method='Nelder-Mead',
			tol=1e-8, options={'initial_simplex':[[0],[1]]})
		epsP = res['x']
		sig[i] = sigma(eps[i],epsP)
		dmg[i] = d(epsP)
		epsPArray[i] = epsP
		epsPOld = epsP

	return sig,dmg,epsPArray

def damagePlasticViscoElastic1D(params, eps, dt):
	"""
	Input: Parameter array, strain array
	Output: stress array
	"""

	# Parameters
	Einf		= params[0]
	C 			= params[1]
	gamma 		= params[2]
	s 			= params[3]
	sig0 		= params[4]
	tau0 		= params[5]
	epsDot0 	= params[6]
	m 			= params[7]
	E			= params[8::2]
	tau			= params[9::2]

	# Exponential damage law 
	def d(epsP):
		return 1-np.exp(-s*abs(epsP))

	# Define stress
	def sigma(eps_i,epsP_i,dt_i,epsP_ve_i):
		return (1-d(epsP_i))*(Einf*(eps_i-epsP_i)+(sum(E_j/(1+dt_i/tau_j)*(eps_i-epsP_ve_ij) for E_j, tau_j, epsP_ve_ij in zip(E, tau, epsP_ve_i))))
		# return ((1-d(epsP_i))*Einf*(eps_i-epsP_i)+(sum(E_j/(1+dt_i/tau_j)*(eps_i-epsP_ve_ij) for E_j, tau_j, epsP_ve_ij in zip(E, tau, epsP_ve_i))))


	# Initialize
	nSteps = len(eps)
	sig = np.zeros(nSteps)
	epsPArray = np.zeros(nSteps)
	epsP_ve_Array = np.zeros((nSteps,len(E)))
	dmg = np.zeros(nSteps)
	# epsPOld = 0
	# epsPveOld = np.zeros(len(E))
	epsPOldepsPveOld=np.zeros(len(E)+1)

	# Integrate
	for i in range(nSteps):
		 
		# Define effective potential
		def f(epsPepsPve):
			epsPOld=epsPOldepsPveOld[0]
			epsPveOld=epsPOldepsPveOld[1:]
			epsP=epsPepsPve[0]
			epsPve=epsPepsPve[1:]
			term1 = 0.5*(1-d(epsP))*Einf*(eps[i]-epsP)*(eps[i]-epsP)
			term10= sum(0.5*(1-d(epsP))*E_j/(1+dt[i]/tau_j)*(eps[i]-epsPve_j)*(eps[i]-epsPve_j) for E_j, tau_j, epsPve_j in zip(E,tau,epsPve))
			# term10= sum(0.5*E_j/(1+dt[i]/tau_j)*(eps[i]-epsPve_j)*(eps[i]-epsPve_j) for E_j, tau_j, epsPve_j in zip(E,tau,epsPve))
			term2 = C/(gamma+1)*pow(np.abs(epsP),gamma+1)
			term3 = sig0 * abs(epsP-epsPOld)
			term4 = tau0/(m+1) * pow(abs(epsP-epsPOld),m+1)/pow(dt[i]*epsDot0,m)
			term40= sum((1-d(epsP))*E_j*tau_j/(2*(1+dt[i]/tau_j))*pow(epsPve_j-epsPveOld_j,2)/dt[i] for E_j,tau_j,epsPve_j,epsPveOld_j in zip(E,tau,epsPve,epsPveOld))
			# term40= sum(E_j*tau_j/(2*(1+dt[i]/tau_j))*pow(epsPve_j-epsPveOld_j,2)/dt[i] for E_j,tau_j,epsPve_j,epsPveOld_j in zip(E,tau,epsPve,epsPveOld))
			return term1 + term10 + term2 + term3 + term4 + term40

		# Define effective potential
		def f_(epsPepsPve):
			epsPOld=epsPOldepsPveOld[0]
			epsPveOld=epsPOldepsPveOld[1:]
			epsP=epsPepsPve[0]
			epsPve=epsPepsPve[1:]
			term1 = 0.5*(1-d(epsP))*Einf*(eps[i]-epsP)*(eps[i]-epsP)
			term10= sum(0.5*(1-d(epsP))*E_j/(1+dt[i]/tau_j)*(eps[i]-epsPve_j)*(eps[i]-epsPve_j) for E_j, tau_j, epsPve_j in zip(E,tau,epsPve))
			# term10= sum(0.5*E_j/(1+dt[i]/tau_j)*(eps[i]-epsPve_j)*(eps[i]-epsPve_j) for E_j, tau_j, epsPve_j in zip(E,tau,epsPve))
			term2 = C/(gamma+1)*pow(np.abs(epsP),gamma+1)
			term3 = sig0 * abs(epsP-epsPOld)
			term4 = tau0/(m+1) * pow(abs(epsP-epsPOld),m+1)/pow(dt[i]*epsDot0,m)
			term40= sum((1-d(epsP))*E_j*tau_j/(2*(1+dt[i]/tau_j))*pow(epsPve_j-epsPveOld_j,2)/dt[i] for E_j,tau_j,epsPve_j,epsPveOld_j in zip(E,tau,epsPve,epsPveOld))
			# term40= sum(E_j*tau_j/(2*(1+dt[i]/tau_j))*pow(epsPve_j-epsPveOld_j,2)/dt[i] for E_j,tau_j,epsPve_j,epsPveOld_j in zip(E,tau,epsPve,epsPveOld))
			print(term1,term10,term2,term3,term4,term40)
			return term1 + term10 + term2 + term3 + term4 + term40

		# Solve for new plastic strain
		initial_simplex = np.zeros((len(E) + 2, len(E) + 1))
		np.fill_diagonal(initial_simplex[1:], 1)
		# res = minimize(f,x0=epsPOldepsPveOld, method='Nelder-Mead',
		# 	tol=1e-8, options={'initial_simplex':initial_simplex})
		x0=epsPOldepsPveOld.copy()
		if i>0:
			x0[0]+=(eps[i]-eps[i-1])
		else:
			x0[0]+=eps[i]
		res = minimize(f,x0, method='Powell',
			tol=1e-8)
		# print(res.nit)
		# print(res.fun)
		# print(f([0.0567619, 0.00894514, 0.0166091]))
		# print(f([0.05354992, -0.01829159]))
		# print(f([0.02336554, 0.01829159]))
		epsPepsPve = res['x']
		epsP=epsPepsPve[0]
		epsPve=epsPepsPve[1:]
		sig[i] = sigma(eps[i],epsP,dt[i],epsPve)
		dmg[i] = d(epsP)
		epsPArray[i] = epsP
		epsP_ve_Array[i] = epsPve
		epsPOldepsPveOld = np.concatenate(([epsP],epsPve))

	return sig,dmg,epsPArray#,epsP_ve_Array

# Parameters
Einf		= 166
C 			= 500
gamma 		= 1.25
s 			= 16
sig0 		= 0
tau0 		= 5
epsDot0 	= 0.007
m 			= 0.25
E			= [447,259,216]
tau			= [6.96,110,1255]

params =[Einf, C, gamma, s, sig0, tau0, epsDot0, m] + E + tau

eps = [0,0.025,0.05,0.075,0.1]
dt = [500,500,500,500,500]

sig,dmg,epsPArray = damagePlasticViscoElastic1D(params, eps, dt)

# Plot the stress vs strain
plt.figure(figsize=(6, 4))
plt.plot(eps, sig, 'b-', label='Stress-Strain')
plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title('Stress vs Strain')
plt.legend()
plt.tight_layout()

# # Plot the damage vs strain
# plt.figure(figsize=(6, 4))
# plt.plot(eps, dmg, 'r-', label='Damage-Strain')
# plt.xlabel('Strain')
# plt.ylabel('Damage')
# plt.title('Damage vs Strain')
# plt.legend()
# plt.tight_layout()

# # Plot the plastic strain vs strain
# plt.figure(figsize=(6, 4))
# plt.plot(eps, epsPArray, 'g-', label='Plastic Strain-Strain')
# plt.xlabel('Strain')
# plt.ylabel('Plastic Strain')
# plt.title('Plastic Strain vs Strain')
# plt.legend()
# plt.tight_layout()

plt.show()
# Plot dynamics of event-driven stochastic processes

import numpy as np
import matplotlib.pyplot as plt
from scipy import linspace

type=0

if type==0:
	result = np.loadtxt('result-sir-net-1') 
elif type==1:
	result = np.loadtxt('result-sir-net-2') 
elif type==2:
	result = np.loadtxt('result-voter-net-N100') 
elif type==3:
	result = np.loadtxt('result-voter-net-N1000') 
elif type==4: # SIR on a metapopularion network
	result = np.loadtxt('result-sir-metapop')
elif type==5: # Lotka-Volterra run 1
	result = np.loadtxt('result-lv-wellmixed-1') 
elif type==6: # Lotka-Volterra run 2
	result = np.loadtxt('result-lv-wellmixed-2') 
	
# First column = time
# Second coumn = fraction of nodes in opinion A

if type==0 or type==1: # SIR
	plt.plot(result[:,0], result[:,1], 'b-', label='S', linewidth=1)
	plt.plot(result[:,0], result[:,2], 'r-', label='I', linewidth=1)
	plt.plot(result[:,0], result[:,3], '-', label='R', linewidth=1, color='brown')
	plt.ylabel('Fraction of individuals', fontsize=20, family='Times New Roman')
	plt.ylim(0, 1)
	plt.yticks(np.arange(0, 1.2, step=0.2))
#	plt.ylim(0, 1)
#	plt.yticks(np.arange(0, 1.2, step=0.2))
elif type==2 or type==3: # voter model
	plt.plot(result[:,0], result[:,1], 'b-', linewidth=1)
	plt.ylabel('Fraction of nodes in opinion A', fontsize=20, family='Times New Roman')
	plt.ylim(0, 1)
	plt.yticks(np.arange(0, 1.2, step=0.2))
elif type==5 or 6: # Lotka-Volterra
	plt.plot(result[:,0], result[:,1], 'b-', label='rabbit', linewidth=1)
	plt.plot(result[:,0], result[:,2], 'r-', label='fox', linewidth=1)
	plt.ylabel('Number of individuals', fontsize=20, family='Times New Roman')
	plt.ylim(0, 2280)
	plt.yticks(np.arange(0, 2500, step=500))

plt.subplots_adjust(bottom=0.16, left=0.2, right=0.95)
plt.xlabel(r'$t$', fontsize=24)
plt.xlim(0, np.nanmax(result[:,0])*1.02) # result[-1,0] = time when the dynamics have terminated
plt.tick_params(labelsize=20)

if type==0:
    plt.title('(a)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')
    plt.legend(loc = 'center right', numpoints = 1, labelspacing=0.25)
elif type==1:
    plt.title('(b)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')
elif type==2:
    plt.title('(a)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')
elif type==3:
    plt.title('(b)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')
elif type==4:
    plt.legend(loc = 'center right', numpoints = 1, labelspacing=0.25)
elif type==5:
    plt.title('(a)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')
    plt.legend(loc = 'upper left', numpoints = 1, labelspacing=0.25)
elif type==6:
    plt.title('(b)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')


# plt.show()
plt.savefig("fig.pdf")
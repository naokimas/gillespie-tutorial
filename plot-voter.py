# Plot dynamics of event-driven stochastic processes

import numpy as np
import matplotlib.pyplot as plt
from scipy import linspace

result = np.loadtxt('result-voter-net') 
# First column = time
# Second coumn = fraction of nodes in opinion A

plt.plot(result[:,0], result[:,1], 'b-', linewidth=1)
plt.subplots_adjust(bottom=0.16, left=0.2, right=0.95)
plt.xlabel(r'$t$', fontsize=24)
plt.ylabel('Fraction of nodes in opinion A', fontsize=20, family='Times New Roman')
plt.xlim(0, np.nanmax(result[:,0])*1.02) # result[-1,0] = time when the dynamics have terminated
plt.ylim(0, 1)
plt.yticks(np.arange(0, 1.2, step=0.2))
plt.tick_params(labelsize=20)

# if data==1:
#    plt.title('(a)', fontsize=28, x=-0.15, y=1.03, family='Times New Roman')
#    plt.legend(loc = 'lower left', numpoints = 1, labelspacing=0.25)

# plt.show()
plt.savefig("fig.pdf")
# Plot change in e (vapor pressure) and esi (saturation vapor pressure wrt ice) for given changes in
# water vapor mixing ratio and temperature 

import numpy as np
import matplotlib.pyplot as plt

savefig = 1
v = 'RHi'
#P = [25000,30000] # Pa
#delta_Q = [4e-6, 1e-5] # kg/kg
#delta_T = [1, 2] 
# Based on case vertcross plots:
#P = [30000, 20000] # Pa
#delta_Q = [1.6e-05, 2.5e-07] # kg/kg
#delta_T = [2, 1]
#color = ['black','gray']
# To demonstrate point best
P = [25000,25000] # Pa
delta_Q = [0.00005,0.00005] # kg/kg
delta_T = [1,2]
color = ['grey','black']


# e:
Q1 = 0.0001 #kg/kg

delta_e = []
for i in range(len(P)):
    #e1 = Q1/(Q1+0.622) * P[i]
    e1 = Q1/0.622 * P[i]
    #e2 = (Q1+delta_Q[i])/((Q1+delta_Q[i])+0.622) * P[i]
    e2 = (Q1+delta_Q[i])/0.622 * P[i]

    delta_e.append(e2 - e1)

# esi:
T1 = np.arange(-50,-20)

plt.figure(figsize=(4,3))
for i in range(len(P)):
    delta_es = []
    for T in T1:
        if v == 'RHw':
            es1 = 6.112*100*np.exp(17.62*T/(243.12+T))
            es2 = 6.112*100*np.exp(17.62*(T+delta_T[i])/(243.12+(T+delta_T[i])))
        elif v == 'RHi':
            es1 = 6.112*100*np.exp(22.46*T/(272.62+T))
            es2 = 6.112*100*np.exp(22.46*(T+delta_T[i])/(272.62+(T+delta_T[i])))

        delta_es.append(es2 - es1)

    plt.plot(T1, delta_es, color=color[i], label='$dT$ %i$\degree$C' %delta_T[i])
    if i == 1:
        plt.plot(T1, np.repeat(delta_e[i],len(T1)),'--', color=color[i])

plt.xlabel('$T_{past} (\degree$C)')
plt.ylabel('$de_{si} / de$ (Pa)')
plt.grid()
plt.legend()
plt.title('$dQ$ 0.05 g/kg')
if savefig:
    plt.savefig('figures/de_desi_RH.png',bbox_inches='tight',dpi=120)
else:
    plt.show(block=False)

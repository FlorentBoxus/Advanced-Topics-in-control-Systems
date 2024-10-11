# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 09:26:22 2023

@author: Florent Boxus
"""

import matplotlib.pyplot as plt
import numpy as np 
import math
from scipy.optimize import fsolve
import sys
##########################################parameters definition
c=1.82
Vr=3
L=0.05## total length
N=400## mesh interval for x-axis
t_step=0.05## time step
ae=6##excitatory
ai=5##inhibitory
re=0.0005##excitatory spread
ri=0.001##inhiboty spread
v=0.08# propagation speed
muP=2.5
alpha1 = 400 
alpha2 = alpha1
S_max=1
############################################adimensionnalisation
L = L/re
v = v/(re*np.sqrt(alpha1*alpha2))

##############################################definition
r=re/ri
alpha = np.sqrt(alpha1/alpha2)
delta_U=L/N
##################################################


def solve_IC(x):
    return (ae -ai)*(1/(1+np.exp(-c*(x-Vr)))) -x + muP


def G(x,u,t):
    space_index1=(x-u)/delta_U
    if space_index1<0:
        space_index1+=N+1
    space_index2=(x+u)/delta_U
    if space_index2>=N+1:
        space_index2-=N+1  
    time_index=(t-u/v)
    if(time_index<0):
        return 2*F_V0
    elif (time_index % t_step)==0:
        return F[int(space_index1),int(time_index/t_step)] + F[int(space_index2),int(time_index/t_step)]
    else:
        
        return F[int(space_index1),int(round(time_index/t_step))] + F[int(space_index2),int(round(time_index/t_step))]
    
    
    
def I(space_index,time_index,type_):
    if type_==1:#excitatory
        eta=1
    else:#inhibitory
        eta=r
    I=0
    for n in range(0,N-1,1):
        I+=(1/delta_U*eta**2)*math.exp(-eta*x[n+1])*(G(x[space_index],x[n+2],t[time_index])-2*G(x[space_index],x[n+1],t[time_index])+G(x[space_index],x[n],t[time_index]))
    I+=(1/eta)*(G(x[space_index],0,t[time_index])+(G(x[space_index],delta_U,t[time_index])-G(x[space_index],0,t[time_index]))/(eta*(delta_U)))
    I-=(1/eta)*math.exp(-eta*L)*(G(x[space_index],L,t[time_index])+(G(x[space_index],L,t[time_index])-G(x[space_index],L-delta_U,t[time_index]))/(eta*(delta_U)))
    return I


def solve_J(index):
    for i in range(0,V.shape[0],1):#i spatial index, index time index
        J[i,index]=(ae/2) * I(i,index,1)- (ai/2)*r * I(i,index,2)
    return


def system_equation(index):
    dphi_dt=np.zeros(V.shape[0])#solve for 1 time step, for all the spatial index
    for i in range(0,V.shape[0],1):
        dphi_dt[i]=-2*phi[i,index]-V[i,index]+J[i,index]+muP
    return dphi_dt



def solve_system_ode():
    for i in range(0,V.shape[1]-1 ,1):#commence à time_step qui est le temps de la perturbation
        solve_F(i)#renvoie les valuers des sigmoides pour tout un pas de temps
        solve_J(i)#renvoie toutes les valeurs de J pour un pas de temps
        dphi_dt=system_equation(i)
        for j in range(0,V.shape[0],1):
            V[j,i+1]=V[j][i]+phi[j,i]*t_step#j is spatial index, i is time index
            phi[j,i+1]=phi[j][i]+dphi_dt[j]*t_step
        print(i)
    return
        
        
def solve_F(index):#prend en argument toutes les valeurs de tension pour un pas de temps
    for i in range(0,V.shape[0],1):
        F[i,index]=S_max/(1+math.exp(-c*(V[i][index]-Vr)))
    return



t=np.arange(0,100+t_step, t_step)
space=L
number_x_steps=int(space/delta_U)#how much mesh grid do we need
x=np.arange(0,space+delta_U,delta_U)#spatial distribution
V0 = fsolve(lambda x: solve_IC(x), 0)#random IC between 0 and 1
V0_eq=V0[0]#when equilibrium
F_V0=S_max/(1+math.exp(-c*(V0_eq-Vr)))
IC=np.random.normal(loc=V0[0], scale=10**-6 * V0[0], size=(N+1, 1))
nbr_t_steps=len(t)
V=np.zeros((len(x),nbr_t_steps))
phi=np.zeros((len(x),nbr_t_steps))
F=np.zeros((len(x),nbr_t_steps))#défini comme variable globale
J=np.zeros((len(x),nbr_t_steps))

for i in range(0,N+1,1):
    V[i,0]=IC[i]#introduction d'une perturbation
solve_system_ode()
np.savetxt('t.txt',t)
np.savetxt('x.txt',x)
np.savetxt('V.txt',V)
plt.pcolormesh(x, t, V, cmap='viridis')  
plt.colorbar(label='V')

# Set labels for x and y axes
plt.xlabel('x')
plt.ylabel('t')

# Add a title if needed
plt.title('Color Plot of V as a Function of x and t')

# Show the plot
plt.show()











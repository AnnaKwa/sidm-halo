from int_functions import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#these params set in [int_functions.py]

#n=100         #number shells
#Mtot=100     #total dimensionless mass, tilde_M= M/M_0
#a=2.257      #order unity constant
#sigma_const=.1  #dimensionless cross section per unit mass, for now will set const but later will be function of r

#conc=4.      #initial NFW halo concentration
#rmin=0.01   #dimensionless radius, tilde_r= r/rs
#rmax=conc*10
#r_arr=np.logspace(np.log10(rmin),np.log10(rmax),num=n,endpoint=True,base=10)

       


########################
#                      #
#  initialize profile  #
#                      #
########################

#initial cross section as function of r
sigma_arr=[]
for i in range(0,n):
    sigma_arr.append(sigma_const)

#initial mass (contained within r_i)
m_arr=[]
for i in range(0,n):
    ri=r_arr[i]
    m_arr.append(-ri/(1+ri) + np.log(1+ri))

#NFW initial density profile
rho_arr=[]
rho_arr.append(3*m_arr[0]/(r_arr[0]**3))
for i in range(1,n):
    rho_arr.append(3*(m_arr[i]-m_arr[i-1])/(r_arr[i]**3-r_arr[i-1]**3) )

#initial pressure profile
p_arr=[]
p_integrand= lambda r: (np.log(1+r)-(r/(1+r)))*(1/(r*(r+1)**2)/r**2) 
p_arr.append(integrate.quad(p_integrand, 0.5*r_arr[0], 100)[0])
for i in range(1,n):
    p_arr.append(integrate.quad(p_integrand, 0.5*(r_arr[i-1]+r_arr[i]), 100)[0])

#initial specific energy profile
u_arr=[]
for i in range(0,n):
    u_arr.append(1.5*p_arr[i]/rho_arr[i])
print rho_arr[0],rho_arr[1]

#initial luminosity profile
L_arr=[]
for i in range(0,n-1):
	L_arr.append(-0.5*(r_arr[i]**2)*np.sqrt((u_arr[i+1]+u_arr[i])/2)*((a*sigma_arr[i]**2 + 2/(p_arr[i]+p_arr[i+1]))**(-1))*2*(u_arr[i+1]-u_arr[i])/(r_arr[i+1]-r_arr[i])  )
#    L_arr.append(-np.sqrt(2./3.)*0.5*r_arr[i]*(np.sqrt(u_arr[i])+np.sqrt(u_arr[i+1]))*(a*sigma_arr[i]**2 + 2/(p_arr[i]+p_arr[i+1]))**(-1) *(u_arr[i+1]-u_arr[i])/(r_arr[i+1]-r_arr[i]) )
#need to estimate nth L point using dL/dr from previous step
L_arr.append(deltaX(L_arr))

print L_arr[0],L_arr[1],L_arr[2]

####################################
#                                  # 
#  plot and save initial profiles  #
#                                  # 
####################################                                 

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_yscale('log',nonposy='clip')
ax.set_xscale('log')
ax.set_xlim([r_arr[1],rmax])
#ax.set_ylim([np.min(p_arr), 2*np.max(p_arr)])

L_abs=[]
init_L_abs=[]
init_u_arr=u_arr
init_rho_arr=rho_arr
init_p_arr=p_arr
for i in range(0,len(r_arr)):
    init_L_abs.append(abs(L_arr[i]))

ax.plot(r_arr, init_u_arr, '--', color='black', label='specific energy init')
ax.plot(r_arr, init_rho_arr, '--', color='blue', label='density init')
ax.plot(r_arr, init_L_abs, '--', color='red', label='luminosity abs. val init')
ax.plot(r_arr, init_p_arr, '--', color='green', label='pressure init')



######################################################################################
#                                                                                    #
#  alternate timesteps in conduction/diffusion and readjusting for hydrostatic eq.   #  
#                                                                                    #
######################################################################################



for numsteps in range(0,6):
    conduction_diffusion(L_arr,m_arr,u_arr,p_arr,rho_arr,r_arr)
    hydroeq_adjust_p(p_arr,rho_arr,u_arr)
    deltaR_arr=tdma(r_arr, p_arr, rho_arr, m_arr,n)
    #deltaR_arr2=linalg_solver(r_arr, p_arr, rho_arr, m_arr,n)

    deltaP_arr,deltaRho_arr=[],[]
    deltaP_arr.append(-5.*p_arr[0]*(r_arr[0]**2)*deltaR_arr[0]/r_arr[0]**3 )
    deltaRho_arr.append(-3.*rho_arr[0]*(r_arr[0]**2)*deltaR_arr[0]/r_arr[0]**3)
    for i in range(0,n):
        if i>0:
            deltaP_arr.append(-5.*p_arr[i]*( (r_arr[i]**2)*deltaR_arr[i] - (r_arr[i-1]**2)*deltaR_arr[i-1] )/(r_arr[i]**3-r_arr[i-1]**3)   )
            deltaRho_arr.append(-3.*rho_arr[i]*( (r_arr[i]**2)*deltaR_arr[i] - (r_arr[i-1]**2)*deltaR_arr[i-1] )/(r_arr[i]**3-r_arr[i-1]**3)   )
    for i in range(0,n):
        p_arr[i]=p_arr[i]+deltaP_arr[i]
        #if p_arr[i]<0.:
        #	p_arr[i]=0.
        rho_arr[i]=rho_arr[i]+deltaRho_arr[i]
        r_arr[i]=r_arr[i]+deltaR_arr[i]
        u_arr[i]=1.5*p_arr[i]/rho_arr[i]
        if p_arr<0.:
        	print 'Error: r['+repr(i)+'], step '+repr(numsteps)+', p<0'
    calcL(L_arr, u_arr, p_arr, r_arr)

for i in range(0,len(r_arr)-1):
	if r_arr[i+1]<r_arr[i]:	
		print 'Error: r['+repr(i+1)+']>r['+repr(i)+']'

for i in range(0,len(r_arr)):
    L_abs.append(abs(L_arr[i]))

#testing purposes
#ax.plot(r_arr, u_arr, '-', color='black', label='specific energy')
ax.plot(r_arr, rho_arr, '-', color='blue', label='density')
ax.plot(r_arr, L_abs, '-', color='red', label='luminosity abs. val')

ax.set_xlabel('dimensionless radius')
ax.set_ylabel('dimensionless quantities')

Lslope=dLdr(L_arr,r_arr)
for elem in Lslope:
    elem=abs(elem)
Lslope.append(deltaX(Lslope))
#ax.plot(r_arr,Lslope, label='dL/dr')
'''
plt.legend(prop={'size':12},
           loc='lower left',
           borderpad=.5,
           ncol=2,
           numpoints=1,
           columnspacing=.25,
           handletextpad=.15)
'''
plt.show()

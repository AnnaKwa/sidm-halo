################################################
#                                              #  
#  define functions for numerical integration  #
#                                              # 
################################################

import numpy as np
from scipy import integrate
n=300        #number shells
a=2.257      #order unity constant
sigma_const=10 #dimensionless cross section per unit mass, for now will set const but later will be function of r
Mtot=1

show_errs=0

#initial cross section as function of r
sigma_arr=[]
for i in range(0,n):
    sigma_arr.append(sigma_const)
#initialize radius arr    
conc=4.      #initial NFW halo concentration
#rmin=0.01   #dimensionless radius, tilde_r= r/rs
rmin=0.05
rmax=conc*5
r_arr=np.logspace(np.log10(rmin),np.log10(rmax),num=n,endpoint=True,base=10)

#returns delta_r*dX/dr from previous shell, add to previous val to get last entry in array
def deltaX(X):
    dXdr=(X[len(X)-1]-X[len(X)-2]) / (r_arr[len(X)-1]-r_arr[len(X)-2])
    deltaX=dXdr*(r_arr[n-1]-r_arr[n-2])
    return deltaX

#calculate time step: want to be small compared to relaxation timescale delta u / u
#returns 0.1*(delta u / u) 
def timestep(uarr,rarr,marr,larr):
    minstep=1e6
    eps_t=1e-3        #fraction delta_u/u
    xmin=0.              #for testing purposes
    for i in range(2,n-3):
    	trelax=-(marr[i+1]-marr[i])/(larr[i+1]-larr[i]) *eps_t*uarr[i]
    	if trelax<minstep and trelax>0.:
    		minstep=trelax

        '''
        if uarr[i]>0.:
            trelax=abs((marr[i+1]-marr[i])/(larr[i+1]-larr[i]))*eps_t*uarr[i]
            if abs(trelax)<minstep:
                minstep=trelax
                xmin=rarr[i]
    	'''
    #print minstep
    return (minstep)

#calculate cooling loss at each r-- will turn this on eventually, for now set to zero
def coolingloss(r):
    return 0.

#calculate luminosity profile
def calcL(Larr,uarr,parr,rarr):
    for i in range(0,n-1):
    	if uarr[i]>0.:
        	#Larr[i]=(-np.sqrt(2./3.)*0.5*(rarr[i]**2)*( np.sqrt(uarr[i])+np.sqrt(uarr[i+1]) )*(a*sigma_arr[i]**2 + 2/(parr[i]+parr[i+1]) )**(-1) *(uarr[i+1]-uarr[i])/(rarr[i+1]-rarr[i]) )
        	Larr[i]=-(rarr[i]**2)*np.sqrt((uarr[i+1]+uarr[i])/3)* ((a*sigma_arr[i]**2 + 2/(parr[i]+parr[i+1]) )**(-1))* 2*(uarr[i+1]-uarr[i])/(rarr[i+1]-rarr[i])
        	#print repr(a*sigma_arr[i]**2 ), repr(2/(parr[i]+parr[i+1]))
        else:
        	#Larr[i]=0.
        	if show_errs==1:
        		print 'Error: u['+repr(i)+']<0'

    Larr[n-1]=(deltaX(Larr)+Larr[n-2])
    

#calculate delta u for each shell and adjust u_i's and p_i's accordingly        
def conduction_diffusion(arrL,arrM,arrU,arrP,rhoarr, rarr):
    delta_u=[]
    delta_p=[]
    delta_t=timestep(arrU,rarr,arrM,arrL)
    for i in range(0,n-1):
        delta_u.append(-delta_t*(arrL[i+1]-arrL[i])/(arrM[i+1]-arrM[i]))
    delta_u.append(delta_u[-1]+deltaX(delta_u))
    
    for i in range(0,n):
        arrU[i]=arrU[i]+delta_u[i]
        if arrU[i]<0.:
            arrU[i]=0.
#    for i in range(0,n-1):
#        delta_p.append((2./3.)*(rhoarr[i]*arrU[i]))
#    delta_p.append(deltaX(arrP))
#    for i in range(0,n):
#        arrP[i]=(2./3)*arrU[i]*rhoarr[i]

def hydroeq_adjust_p(arrP,arrRho,arrU):
    for i in range(0, n):
        arrP[i]=(2./3)*arrU[i]*arrRho[i]

def dLdr(arrL,arrR):
    dLdr_arr=[]
    for i in range(0,n-1):
        dLdr_arr.append((arrL[i+1]-arrL[i])/(arrR[i+1]-arrR[i]))
    return dLdr_arr

def create_tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def linalg_solver(r,p,rho,M,n):
	#define tridiagonal matrix A	
    a,b,c,d=[],[],[],[]
    for i in range(1,n-1):
    	aa=-10.*(r[i]**2)*(r[i-1]**2)*p[i]/(r[i]**3-r[i-1]**3) + 1.5*M[i]*(r[i-1]**2)*(r[i+1]-r[i-1])*rho[i]/(r[i]**3-r[i-1]**3) - 0.5*M[i]*(rho[i]+rho[i+1])
    	a.append(aa)
    aa=-10*(r[n-1]**2)*(r[n-2]**2)*p[n-1]/(r[n-1]**3-r[n-2]**3) + 1.5*M[n-1]*(r[n-2]**2)*( (r[n-1]+deltaX(r))-r[n-2])*rho[n-1]/(r[n-1]**3-r[n-2]**3) - 0.5*M[n-1]*(rho[n-1]+ rho[n-1]+deltaX(rho))
    a.append(aa)
    bb=10.*(r[0]**4)*( (p[1]/(r[1]**3-r[0]**3)) + (p[0]/(r[0]**3)) ) + 4.*r[0]*(p[1]-p[0]) + 1.5*(r[0]**2)*(r[1])*M[0]*( rho[1]/(r[1]**3-r[0]**3) - rho[0]/r[0]**3  )
    b.append(bb)
    for i in range(1,n-1):
        #b.append( 10.*(r[i]**4)*(p[i+1]/(r[i+1]**3-r[i]**3) - p[i]/(r[i]**3-r[i-1]**3) ) + 4.*r[i]*(p[i+1]-p[i]) + 1.5*(r[i]**2)*(r[i+1]-r[i-1])*M[i]*( rho[i+1]/(r[i+1]**3-r[i]**3) - rho[i]/(r[i]**3-r[i-1]**3) )  )
        bb=10.*(r[i]**4)*(p[i+1]/(r[i+1]**3-r[i]**3) + p[i]/(r[i]**3-r[i-1]**3) ) + 4.*r[i]*(p[i+1]-p[i]) + 1.5*M[i]*( rho[i]*(r[i]**2)*(r[i+1] - r[i-1])/(r[i]**3-r[i-1]**3) + rho[i+1]*(r[i]**2)*(r[i+1] - r[i-1])/(r[i]**3-r[i-1]**3) )
        b.append(bb)
      
    bb=10.*(r[n-1]**4)*( (p[n-1]+deltaX(p))/( (r[n-1]+deltaX(r))**3-r[n-2]**3) + p[n-1]/(r[n-1]**3-r[n-2]**3) ) + 4.*r[n-1]*(p[n-1]+deltaX(p)-p[n-1]) + 1.5*(r[n-1]**2)*M[n-1]*(r[n-2]-(r[n-1]+deltaX(r)) )/(r[n-1]**3-r[n-2]**3) + (rho[n-1]+deltaX(rho))*(r[n-1]+deltaX(r) - r[n-2])/((r[n-1]+deltaX(r))**3-r[n-1]**3)
    b.append(bb)
    cc=-10.*(r[0]**2)*(r[1]**2)*p[1]/(r[1]**3-r[0]**3)-1.5*M[0]*(r[1]**2)*(r[1]-r[0])*rho[1]/(r[1]**3-r[0]**3)+0.5*M[0]*(rho[1]+rho[0])	
    c.append(cc)
    #c.append(-10.*(r[0]**2)*(r[1]**2)*p[1]/(r[1]**3-r[0]**3) - 1.5*M[0]*(r[1]**2)*(r[1]-r[0])*rho[1]/(r[1]**3-r[0]**3) + 0.5*M[0]*(rho[1]+rho[0])  )
    for i in range(1,n-1):
        cc=-10.*(r[i]**2)*(r[i+1]**2)*p[i+1]/(r[i+1]**3-r[i]**3) - 1.5*M[i]*(r[i+1]**2)*(r[i+1]-r[i-1])*rho[i+1]/(r[i+1]**3-r[i]**3) + 0.5*M[i]*(rho[i+1]+rho[i])
        c.append(cc)
    
    A=create_tridiag(a,b,c)

    dd=2.*(r[0]**2)*(p[1]-p[0]) + 0.5*M[0]*r[1]*(rho[1]+rho[0])
    d.append(-dd)
    for i in range(1,n-1):
        dd=2.*(r[i]**2)*(p[i+1]-p[i]) + 0.5*M[i]*(r[i+1]-r[i-1])*(rho[i+1]+rho[i])
        d.append(-dd)
	dd=(2.*(r[n-1]**2)*( (p[n-1]+deltaX(p)) - p[n-1] ) + 0.5*M[n-1]*(r[n-1]+deltaX(r) - r[n-2])*(rho[n-1]+deltaX(rho) + rho[n-1]) )   	   	
    d.append(-dd)

    del(aa,bb,cc,dd)

    deltar_arr=np.linalg.solve(A,d)
    return(deltar_arr)



def tdma(r,p,rho,M,n):
    a,b,c,d=[],[],[],[]
    x=[]
    a.append(0.)
    for i in range(1,n-1):
    	aa=-10.*(r[i]**2)*(r[i-1]**2)*p[i]/(r[i]**3-r[i-1]**3) + 1.5*M[i]*(r[i-1]**2)*(r[i+1]-r[i-1])*rho[i]/(r[i]**3-r[i-1]**3) - 0.5*M[i]*(rho[i]+rho[i+1])
    	a.append(aa)
    aa=-10*(r[n-1]**2)*(r[n-2]**2)*p[n-1]/(r[n-1]**3-r[n-2]**3) + 1.5*M[n-1]*(r[n-2]**2)*( (r[n-1]+deltaX(r))-r[n-2])*rho[n-1]/(r[n-1]**3-r[n-2]**3) - 0.5*M[n-1]*(rho[n-1]+ rho[n-1]+deltaX(rho)) 
    a.append(aa)
    #a.append(-10*(r[n-1]**2)*(r[n-2]**2)*p[n-1]/(r[n-1]**3-r[n-2]**3) + 1.5*M[n-1]*(r[n-2]**2)*( (r[n-1]+deltaX(r))-r[n-2])*rho[n-1]/(r[n-1]**3-r[n-2]**3) - 0.5*M[n-1]*(rho[n-1]+ rho[n-1]+deltaX(rho)) )

    #b.append(10.*(r[0]**4)*( (p[1]/(r[1]**3-r[0]**3)) - (p[0]/(r[0]**3)) ) + 4.*r[0]*(p[1]-p[0]) + 1.5*(r[0]**2)*(r[1])*M[0]*( rho[1]/(r[1]**3-r[0]**3) - rho[0]/r[0]**3  ) ) 
    bb=10.*(r[0]**4)*( (p[1]/(r[1]**3-r[0]**3)) + (p[0]/(r[0]**3)) ) + 4.*r[0]*(p[1]-p[0]) + 1.5*(r[0]**2)*(r[1])*M[0]*( rho[1]/(r[1]**3-r[0]**3) - rho[0]/r[0]**3  )
    b.append(bb)
    for i in range(1,n-1):
        #b.append( 10.*(r[i]**4)*(p[i+1]/(r[i+1]**3-r[i]**3) - p[i]/(r[i]**3-r[i-1]**3) ) + 4.*r[i]*(p[i+1]-p[i]) + 1.5*(r[i]**2)*(r[i+1]-r[i-1])*M[i]*( rho[i+1]/(r[i+1]**3-r[i]**3) - rho[i]/(r[i]**3-r[i-1]**3) )  )
        bb=10.*(r[i]**4)*(p[i+1]/(r[i+1]**3-r[i]**3) + p[i]/(r[i]**3-r[i-1]**3) ) + 4.*r[i]*(p[i+1]-p[i]) + 1.5*M[i]*( rho[i]*(r[i]**2)*(r[i+1] - r[i-1])/(r[i]**3-r[i-1]**3) + rho[i+1]*(r[i]**2)*(r[i+1] - r[i-1])/(r[i]**3-r[i-1]**3) )
        b.append(bb)
      
    bb=10.*(r[n-1]**4)*( (p[n-1]+deltaX(p))/( (r[n-1]+deltaX(r))**3-r[n-2]**3) + p[n-1]/(r[n-1]**3-r[n-2]**3) ) + 4.*r[n-1]*(p[n-1]+deltaX(p)-p[n-1]) + 1.5*(r[n-1]**2)*M[n-1]*(r[n-2]-(r[n-1]+deltaX(r)) )/(r[n-1]**3-r[n-2]**3) + (rho[n-1]+deltaX(rho))*(r[n-1]+deltaX(r) - r[n-2])/((r[n-1]+deltaX(r))**3-r[n-1]**3)
    b.append(bb)
    cc=-10.*(r[0]**2)*(r[1]**2)*p[1]/(r[1]**3-r[0]**3)-1.5*M[0]*(r[1]**2)*(r[1]-r[0])*rho[1]/(r[1]**3-r[0]**3)+0.5*M[0]*(rho[1]+rho[0])	

    c.append(cc)
    #c.append(-10.*(r[0]**2)*(r[1]**2)*p[1]/(r[1]**3-r[0]**3) - 1.5*M[0]*(r[1]**2)*(r[1]-r[0])*rho[1]/(r[1]**3-r[0]**3) + 0.5*M[0]*(rho[1]+rho[0])  )
    for i in range(1,n-1):
        cc=-10.*(r[i]**2)*(r[i+1]**2)*p[i+1]/(r[i+1]**3-r[i]**3) - 1.5*M[i]*(r[i+1]**2)*(r[i+1]-r[i-1])*rho[i+1]/(r[i+1]**3-r[i]**3) + 0.5*M[i]*(rho[i+1]+rho[i])
        c.append(cc)
    c.append(0.)
    del(cc)

    dd=2.*(r[0]**2)*(p[1]-p[0]) + 0.5*M[0]*r[1]*(rho[1]+rho[0])
    d.append(dd)
    for i in range(1,n-1):
        dd=2.*(r[i]**2)*(p[i+1]-p[i]) + 0.5*M[i]*(r[i+1]-r[i-1])*(rho[i+1]+rho[i])
        d.append(dd)
	dd=(2.*(r[n-1]**2)*( (p[n-1]+deltaX(p)) - p[n-1] ) + 0.5*M[n-1]*(r[n-1]+deltaX(r) - r[n-2])*(rho[n-1]+deltaX(rho) + rho[n-1]) )   	   	
    d.append(dd)

    #check if algorithm is stable: requires |b_i|>|a_i|+|c_i| for all i
    for j in range(0,len(a)):
    	if (a[j]+c[j])>b[j]:
    		if show_errs==1:
        		print 'Error: TDMA solver unstable for this matrix for i='+repr(i)

    nf = len(a)     # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
    for it in xrange(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]

    xc = ac
    xc[-1] = dc[-1]/bc[-1]

    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    del bc, cc, dc  # delete variables from memory
    return xc

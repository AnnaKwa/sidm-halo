import numpy as np


def tdma(r,p,rho,M,n):
    a,b,c,d=[],[],[],[]
    a.append(0.)
    for i in range(1,n-1):
        a.append(-10.*(r[i]**2)*(r[i-1]**2)*p[i]/(r[i]**3-r[i-1]**3) + 1.5*M[i]*(r[i-1]**2)*(r[i+1]-r[i-1])*rho[i]/(r[i]**3-r[i-1]**3) - 0.5*M[i]*(rho[i]+rho[i+1]) )
    a.append(-10*(r[n-1]**2)*(r[n-2]**2)*p[n-1]/(r[n-1]**3-r[n-2]**3) + 1.5*M[n-1]*(r[n-2]**2)*( (r[n-1]+deltaX(r,r))-r[n-2])*rho[n-1]/(r[n-1]**3-r[n-2]**3) - 0.5*M[n-1]*(rho[n-1]+ rho[n-1]+deltaX(rho,r)) )

    b.append(10.*(r[0]**4)*( (p[1]/(r[1]**3-r[0]**3)) - (p[0]/(r[0]**3)) ) + 4.*r[0]*(p[1]-p[0]) + 1.5*(r[0]**2)*(r[1])*M[0]*( rho[1]/(r[1]**3-r[0]**3) - rho[0]/r[0]**3  ) ) 
    for i in range(1,n-1):
        b.append( 10.*(r[i]**4)*(p[i+1]/(r[i+1]**3-r[i]**3) - p[i]/(r[i]**3-r[i-1]**3) ) + 4.*r[i]*(p[i+1]-p[i]) + 1.5*(r[i]**2)*(r[i+1]-r[i-1])*M[i]*( rho[i+1]/(r[i+1]**3-r[i]**3) - rho[i]/(r[i]**3-r[i-1]**3) )  )
    b.append( 10.*(r[n-1]**4)*( (p[n-1]+deltaX(p,r))/( (r[n-1]+deltaX(r,r))-r[n-2]) + p[n-1]/(r[n-1]**3-r[n-2]**3) ) + 4.*r[n-1]*(p[n-1]+deltaX(p,r)-p[n-1]) + 1.5*(r[n-1]**2)*M[n-1]*( (rho[n-1]+deltaX(rho,r))/((r[n-1]+deltaX(r,r))**3-r[n-1]**3 ) - rho[n-1]/(r[n-1]**3 - r[n-2]**3) )    )

    c.append(-10.*(r[0]**2)*(r[1]**2)*p[1]/(r[1]**3-r[0]**3) - 1.5*M[0]*(r[1]**2)*(r[1]-r[0])*rho[1]/(r[1]**3-r[0]**3) + 0.5*M[0]*(rho[1]+rho[0])  )
    for i in range(1,n-1):
        c.append(-10.*(r[i]**2)*(r[i+1]**2)*p[i+1]/(r[i+1]**3-r[i]**3) - 1.5*M[i]*(r[i+1]**2)*(r[i+1]-r[i-1])*rho[i+1]/(r[i+1]**3-r[i]**3) + 0.5*M[i]*(rho[i+1]-rho[i])   )
    c.append(0.)

    d.append(2.*(r[0]**2)*(p[1]-p[0]) + 0.5*M[0]*r[1]*(rho[1]+rho[0])   )
    for i in range(1,n-1):
        d.append(2.*(r[i]**2)*(p[i+1]-p[i]) + 0.5*M[i]*(r[i+1]-r[i-1])*(rho[i+1]+rho[i])    )
    d.append(2.*(r[n-1]**2)*( (p[n-1]+deltaX(p,r)) - p[n-1] ) + 0.5*M[n-1]*(r[n-1]+deltaX(r,r) - r[n-2])*(rho[n-1]+deltaX(rho,r) + rho[n-1])  )

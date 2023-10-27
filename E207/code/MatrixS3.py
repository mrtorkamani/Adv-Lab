#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 19:45:49 2023

@author: mr
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:03:25 2023

@author: mr
"""
import sympy as sp
import numpy as np

def poly2(x, a, b, c):
    return a*x**2 + b*x + c


kq = sp.symbols('k')
be, al, ga = sp.symbols('be'), sp.symbols('al'), sp.symbols('ga')

def quadrupole(k, L, dim="z"):
    if dim == "x":
        k = -k
    out = np.array( [ [1, L], [k * L, 1] ] )
    return out

def drift(L):
    return np.array( [ [1, L], [0, 1] ] )

def kickangle(mrad):
    return np.array([[0], [mrad]]).T


#print(MS2z)

Md3 = drift(0)
MQ3x = quadrupole(kq, 0.074, "x")
MQ3z = quadrupole(kq, 0.074, "z")
Md4 = drift(.03+.0785+.316+.02)

tmpx = np.matmul(MQ3x, Md3)
MS3x = np.matmul(Md4, tmpx)

tmpz = np.matmul(MQ3z, Md3)
MS3z = np.matmul(Md4, tmpz)


#print(" ")
print(MS3x)

term1x = sp.expand(MS3x[0,0]*MS3x[0,0]) * be  # eps_beta0
term2x = sp.expand(-2*MS3x[0,0]*MS3x[0,1]) * al # eps_alpha0
term3x = sp.expand(MS3x[0,1]**2) * ga        # eps_gamma0    

term1z = sp.expand(MS3z[0,0]*MS3z[0,0]) * be  # eps_beta0
term2z = sp.expand(-2*MS3z[0,0]*MS3z[0,1]) * al # eps_alpha0
term3z = sp.expand(MS3z[0,1]**2) * ga        # eps_gamma0    


collectedx = sp.collect(sp.expand(term1x + term2x + term3x),kq)
collectedz = sp.collect(sp.expand(term1z + term2z + term3z),kq)

#print(collected)


from sympy.solvers import solve
from sympy import Symbol

pop2x = np.array([1.15655339e-02, -7.42838741e-01,  1.52161729e+01])
pop2z = np.array([7.00093697e-03, 4.89874898e-01, 9.21505063e+00])
perr2x= np.array([0.00098806, 0.02047723, 0.77856152])
perr2z= np.array([0.00026697, 0.00553285, 0.21036349])

pxmax2 = pop2x+perr2x
pxmin2 = pop2x-perr2x
pzmax2 = pop2z+perr2z
pzmin2 = pop2z-perr2z

dx2 = pop2x[0]*kq**2 + pop2x[1]*kq + pop2x[2]
dxmax2 = pxmax2[0]*kq**2 + pxmax2[1]*kq + pxmax2[2]
dxmin2 = pxmin2[0]*kq**2 + pxmin2[1]*kq + pxmin2[2]

dz2 = pop2z[0]*kq**2 + pop2z[1]*kq + pop2z[2]
dzmax2 = pzmax2[0]*kq**2 + pzmax2[1]*kq + pzmax2[2]
dzmin2 = pzmin2[0]*kq**2 + pzmin2[1]*kq + pzmin2[2]

#print(pxmax2, perr2x, pop2x, pxmin2)

solx2 = solve(collectedx-dx2, [al, be, ga], dict=True)
solxmax2 = solve(collectedx-dxmax2, [al, be, ga], dict=True)
solxmin2 = solve(collectedx-dxmin2, [al, be, ga], dict=True)


erralx = (solxmax2[0][al] - solx2[0][al])
errbex = (solxmax2[0][be] - solx2[0][be])
errgax = (solxmax2[0][ga] - solx2[0][ga])
#solxmax2
#solxmin2

print('Screen & Direction & $\epsilon \\alpha_0$$\mathrm{[mm^2m^{-1}]}$ & $\epsilon \\beta_0$$\mathrm{[mm^2]}$ &$\epsilon \\gamma_0$$\mathrm{[mm^2m^{-2}]}$ \\\ ')
print('S2 & x & $%f \pm %f$ & $%f \pm %f$ & $%f \pm %f$ \\\ ' %(solx2[0][al],erralx,solx2[0][be],errbex,solx2[0][ga],errgax))

eps2x = float(solx2[0][be])*float(solx2[0][ga])-float(solx2[0][al])**2
eps2xmax = solxmax2[0][be]*solxmax2[0][ga]-solxmax2[0][al]**2
eps2xmin = solxmin2[0][be]*solxmin2[0][ga]-solxmin2[0][al]**2
epsx2 = np.sqrt(float(eps2x))
epsxmax2 = np.sqrt(float(eps2xmax))
epsxmin2 = np.sqrt(float(eps2xmin))

epsxerr2 = (np.abs(epsx2-epsxmin2) + np.abs(epsx2-epsxmax2))/2

epsxerr_sys2 = np.sqrt(epsxerr2**2 + 0.107**2)

#print('Beam emittance in x:', epsx2, epsxerr2, epsxerr_sys2)#, epsxmax, epsxmin)

solz2 = solve(collectedz-dz2, [al, be, ga], dict=True)
solzmax2 = solve(collectedz-dzmax2, [al, be, ga], dict=True)
solzmin2 = solve(collectedz-dzmin2, [al, be, ga], dict=True)

erralz = (solzmax2[0][al] - solz2[0][al])
errbez = (solzmax2[0][be] - solz2[0][be])
errgaz = (solzmax2[0][ga] - solz2[0][ga])

#print(solz)
#print("TWISS2",solz2, solzmax2, solzmin2)
print('S2 & z & $%f \pm %f$ & $%f \pm %f$ & $%f \pm %f$ \\\ ' %(solz2[0][al],erralz,solz2[0][be],errbez,solz2[0][ga],errgaz))

eps2z = solz2[0][be]*solz2[0][ga]-solz2[0][al]**2
eps2zmax = solzmax2[0][be]*solzmax2[0][ga]-solzmax2[0][al]**2
eps2zmin = solzmin2[0][be]*solzmin2[0][ga]-solzmin2[0][al]**2
epsz2 = np.sqrt(float(eps2z))
epszmax2 = np.sqrt(float(eps2zmax))
epszmin2 = np.sqrt(float(eps2zmin))

epszerr2 = (np.abs(epsz2-epszmin2) + np.abs(epsz2-epszmax2))/2

epszerr_sys2 = np.sqrt(epszerr2**2 + 0.021**2)

print('Beam emittance in x: $%f \pm %f$' %(epsx2, epsxerr2))#, epsxmax, epsxmin)
print('Beam emittance in z: $%f \pm %f$' %(epsz2, epszerr2))#, epszmax, epszmin)






kw_x = - pop2x[1]/(2*pop2x[0])
kw_x_err = np.sqrt((perr2x[1]/(2*pop2x[0]))**2 + (pop2x[1]*perr2x[0]/(2*pop2x[0]**2))**2 )
print("Waist k_x:$",kw_x,'\pm', kw_x_err,'$')

sigw_x = poly2(kw_x, *pop2x)
sigw_x_err = np.sqrt( 0*((2*pop2x[0]*kw_x+pop2x[1])*kw_x_err)**2 + (kw_x**2*perr2x[0])**2 + (kw_x*perr2x[1])**2 + perr2x[2]**2 )
print("Waist sig_x^2:$",sigw_x,'\pm', sigw_x_err,'$')


MQ2x = quadrupole(kw_x, 0.074, "x")


tmpx = np.matmul(MQ2x, Md3)
MS2x = np.matmul(Md4, tmpx)

#print(MS2)

tmp  = sigw_x * ( MS2x[1,0]*MS2x[1,0]*float(solx2[0][be]) - 2*MS2x[1,0]*MS2x[1,1]*float(solx2[0][al]) + MS2x[1,1]**2 * float(solx2[0][ga]) )
epsx = np.sqrt(float(tmp))
epsx_err = sigw_x_err / np.sqrt(sigw_x)
print("eps_x:$", epsx,'\pm', epsx_err,'$')



kw_z = - pop2z[1]/(2*pop2z[0])
kw_z_err = np.sqrt((perr2z[1]/(2*pop2z[0]))**2 + (pop2z[1]*perr2z[0]/(2*pop2z[0]**2))**2 )
print("Waist k_z:$",kw_z,'\pm', kw_z_err,'$')

sigw_z = poly2(kw_z, *pop2z)
sigw_z_err = np.sqrt( 0*((2*pop2z[0]*kw_z+pop2z[1])*kw_z_err)**2 + (kw_z**2*perr2z[0])**2 + (kw_z*perr2z[1])**2 + perr2z[2]**2 )
print("Waist sig_z^2$",sigw_z,'\pm', sigw_z_err,'$')


MQ2z = quadrupole(kw_z, 0.074, "x")

tmpz = np.matmul(MQ2z, Md3)
MS2z = np.matmul(Md4, tmpz)

tmp  = sigw_z * ( MS2z[1,0]*MS2z[1,0]*solz2[0][be] - 2*MS2z[1,0]*MS2z[1,1]*solz2[0][al] + MS2z[1,1]**2 * solz2[0][ga] )
epsz = np.sqrt(float(tmp))
epsz_err = sigw_z_err / np.sqrt(sigw_z)
print("eps_z:$", epsz,'\pm', epsz_err,'$')





















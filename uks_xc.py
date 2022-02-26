from __future__ import division
import numpy as np

def lda_functional(x_name, c_name, rho):
    #exchange-correlation functionals

    epsilon = 1e-50
    rho_sum = rho[0] + rho[1] + epsilon

    ex = ec = vx_a = vx_b = vc_a = vc_b = 0

    if x_name == 'LDA':
        #spin polarized Local Spin Density Approximation

        k = pow(3/np.pi,1/3)
        alpha = 2/3

        #exchange energy
        ex = -9/8 * alpha * k * pow(2,1/3) * (pow(rho[0],4/3) + pow(rho[1],4/3)) / rho_sum

        #exchange potential
        vx_a = ex - 0.375 * alpha * k * pow(2,1/3) * (pow(rho[0],4/3) + 4*pow(rho[0],1/3)*rho[1] - 3*pow(rho[1],4/3)) / rho_sum
        vx_b = ex - 0.375 * alpha * k * pow(2,1/3) * (pow(rho[1],4/3) + 4*pow(rho[1],1/3)*rho[0] - 3*pow(rho[0],4/3)) / rho_sum

    if c_name == 'VWN_RPA':
        def evaluate_vwn_rpa_e(t, x):
            #evalute the expression for VWN RPA correlation energy

            ec  =  t[3] * np.log( x*x/X(x,t[1],t[2]))
            ec +=  t[3] * 2*t[1]/Q(t[1],t[2])* np.arctan(Q(t[1],t[2]) / (2*x + t[1]))
            ec += -t[3] * t[1]*t[0]/X(t[0],t[1],t[2])*np.log(pow(x - t[0],2) / X(x,t[1],t[2]))
            ec += -t[3] * t[1]*t[0]/X(t[0],t[1],t[2])*2*(t[1]+2*t[0])/Q(t[1],t[2])*np.arctan(Q(t[1],t[2]) / (2*x + t[1]))

            return ec

        def evaluate_vwn_rpa_v(t, p, x):
            #evaluate derivative of VWN RPA correlation energy

            vc  = t[3]/3 * (-1 + x /(x*x + t[1] * x + t[2]) *(2*x + t[1])/2 )  / rho_sum 
            vc += p * x/ pow(2*x + t[1],2) *2  / rho_sum  /( pow(Q(t[1],t[2]),2) / pow(2*x + t[1],2) + 1) 
            vc += t[0]*t[1]/(t[0]*t[0]+t[0]*t[1]+t[2]) *t[3]*x/rho_sum/3*(1 - (x-t[0])/(x*x+t[1]*x+t[2])*(x+t[1]/2)) /(x-t[0]) 

            return vc

        #square root Seitz radius
        x = pow(3/(np.pi * rho_sum * 4),1/6)

        #helper lambdas
        Q = lambda b,c:np.sqrt(4*c - b*b)
        X = lambda x,b,c:x*x + b*x + c

        #paramagnetic component
        p = [-0.409286, 13.0720,  42.7198, 0.0310907]
        ecp = evaluate_vwn_rpa_e(p, x)

        #ferromagnetic component
        f = [-0.743294, 20.1231, 101.578, 0.01554535]
        ecf = evaluate_vwn_rpa_e(f, x)

        #spin polarization
        zeta = (rho[0] - rho[1])/rho_sum
        fz = (pow(1+zeta,4/3) + pow(1-zeta,4/3) - 2)/(2 * pow(2,1/3) - 2)

        ec = ecp - fz * (ecp - ecf)

        #paramagnetic derivative with respect to x
        pa = (p[1]*p[3] - p[3]*p[1]*p[0]*(p[1]+2*p[0])/(p[0]*p[0]+p[0]*p[1]+p[2]))/3
        vecp = evaluate_vwn_rpa_v(p, pa, x)

        #ferromagnetic derivative with respect to x
        fa = (f[1]*f[3] - f[3]*f[1]*f[0]*(f[1]+2*f[0])/(f[0]*f[0]+f[0]*f[1]+f[2]))/3
        vecf = evaluate_vwn_rpa_v(f, fa, x)
        
        #derivatives of fz with respect to zeta
        dfz = 4/3*pow(2/rho_sum,1/3)*(pow(rho[0],1/3) - pow(rho[1],1/3)) / (2 * pow(2,1/3) - 2)

        #derivative of zeta with respect to rho alpha
        dza = 2*rho[1]/pow(rho_sum,2)
        vc_a = ec + rho_sum * (vecp*(1-fz) - dfz*dza *(ecp - ecf) + vecf*fz)

        #derivatives with respect to rho beta
        dzb = 4/3*rho[0]/pow(rho_sum,2)*pow(2/rho_sum,1/3)*(pow(rho[1],1/3) - pow(rho[0],1/3))/ ( pow(2,1/3) - 1) 
        vc_b = ec + rho_sum * (vecp*(1-fz) - ecp*dzb  + vecf * fz + ecf*dzb) 

    return ec+ex , (vc_a+vx_a, vc_b+vx_b)

def evaluate_gto(gto, p):
    #compute the value of gaussian density at (x,y,z)

    A = (p - gto.atom.center) ; L = np.prod(A**gto.momenta, axis=1).reshape(-1,1)

    phi = np.sum(L*gto.normals*gto.coefficients*np.exp(-gto.exponents*np.sum(A*A, axis=1).reshape(-1,1)), axis=1)

    return phi.reshape(-1,1)

def evaluate_atomic_orbital(basis, p):
    #evaluate the GTO of the atomic orbitals of the molecule

    ao = []
    for i in basis:
        ao.append(evaluate_gto(i, p))

    return np.hstack(ao)

def evaluate_rho_lda(d, ao, weights):
    #evaluate the density over grid shells

    d = d + np.transpose(d, axes=(0,2,1))
    ao_density = np.einsum('pr,xrq->xpq', ao, d, optimize=True)
    ao_density = np.einsum('pi,xpi->xp', ao, ao_density, optimize=True)

    #set small values to zeros
    ao_density[abs(ao_density) < 1.0e-15] = 0

    return ao_density

def evaluate_vxc(vrho, ao, weights):
    #construct exchange-correlation matrix

    weighted_ao = np.einsum('pi,xp->xpi', ao, 0.5*weights*vrho, optimize=True)
    vxc = np.einsum('rp,xrq->xpq', ao, weighted_ao, optimize=True)
    vxc = vxc + np.transpose(vxc, axes=(0,2,1))

    return vxc

def evaluate_exc(exc, rho, weights):
    #evaluate exchange-correlation energy

    return np.einsum('xp,p->x', rho*weights, exc, optimize=True)

from __future__ import division
import numpy as np

from uks_grid import GRID
from uks_xc import evaluate_atomic_orbital, evaluate_rho_lda, lda_functional, evaluate_vxc, evaluate_exc
from uks_util import dii_subspace, out

class pKS(object):
    #polarized Kohn-Sham - performs a KS computation
    '''
    properties - inherits molecules
               - xc        [string]
               - cycles    [integer]
               - tol       [float]
               - mesh      [string]

               - diis      [boolean]
               - diis_size [integer]
    '''
    def __init__(self, mol, xc, cycles=50, tol=1e-6, mesh='close', diis=True, diis_size=6):

        self.mol = mol
        self.x, self.c = xc.split('&')
        self.cycles = cycles
        self.tol=tol
        self.mesh = mesh
        self.mom = self.get_orbital_occupation

        #do consistency checks 
        if ((self.mol.spin + np.sum(self.mol.nele)) % 2) != 0:
            print('spin and number of electrons incompatible')
        self.open = (self.mol.spin != 0)

        self.DIIS = diis ; self.diis_size = diis_size

        #default is to converge to restricted solution
        self.closed_shell_behavior = 'r'

        #output details of molecule
        out(self.mol.silent, [self.mol, self.DIIS, self.diis_size, xc, self.mesh, self.open], 'initial')
        out(self.mol.silent, [cycles, tol], 'cycle')

        self._cache = {}
        
    def get(self, key):
        #retrieve computed values

        return self._cache[key]

    def closed_shell(self, mo_coeff):
        #break symmetry for unrestriced closed-shell
        from math import cos, sin

        #rotate homo and lumo by 45 degrees
        homo = self.mol.nele[0]
        theta = 0.25 * np.pi

        #rotation matrix
        rotate = np.array([[cos(theta), -sin(theta)],
                           [sin(theta),  cos(theta)]])

        #rotate homo and lumo
        c = mo_coeff[0][:, homo-1:homo+1]
        mo_coeff[0][:,homo-1:homo+1] = np.einsum('pr,qr->qp',  rotate, c, optimize=True)

        return mo_coeff

    def get_density_matrix(self, mo_coeff, mo_occ):
        #construct the one electron density matrix

        da =  np.einsum('pr,qr->pq', mo_coeff[0]*mo_occ[0], mo_coeff[0], optimize=True)
        db  = np.einsum('pr,qr->pq', mo_coeff[1]*mo_occ[1], mo_coeff[1], optimize=True)

        return np.array([da, db])

    def get_spin_statistics(self, mo_coeff, mo_occ, s):
        #compute the actual spin squared and multiplicity

        alpha, beta = (mo_coeff[0][:, mo_occ[0]>0], mo_coeff[1][:,mo_occ[1]>0])
        occupations = (alpha.shape[1], beta.shape[1])

        s = np.einsum('rp,rs,sq->pq', alpha, s, beta, optimize=True)

        #get spin components [xy] and [z]
        spin = []
        spin.append(sum(occupations) * 0.5 - np.einsum('ij,ij->', s, s))
        spin.append(pow(occupations[1] - occupations[0], 2) * 0.25)
        ss = sum(spin)

        return ss, 2 * np.sqrt(ss + 0.25)


    def get_orbital_occupation(self, mo_energy, *_):
        #determine occupation numbers of orbitals

        #sort eigenvalues
        e_idx  = (np.argsort(mo_energy[0]), np.argsort(mo_energy[1]))
        e_sort = (mo_energy[0][e_idx[0]], mo_energy[1][e_idx[1]])

        #occupied orbitals - set to 1 the lowest occupied
        mo_occ = np.zeros_like(mo_energy)
        mo_occ[0, e_idx[0][:self.mol.nele[0]]] = 1 ;  mo_occ[1, e_idx[1][:self.mol.nele[1]]] = 1

        return mo_occ

    def maximum_overlap_method(self, occorb, setocc):
        #The Maximum Overlap Method

        imposed_occupation = (occorb[0][:, setocc[0]>0] , occorb[1][:, setocc[1]>0])

        def get_orbital_occupation(_ , mo_coeff, s):
            #Use the MoM to generate next generation of occupation numbers
 
            mo_occ = np.zeros_like(setocc)

            for spin in [0, 1]:

                #compute the overlap between old and new coefficients
                occupation_number = int(np.sum(setocc[spin]))
                mom_s = np.einsum('rp,rs,sq->pq', imposed_occupation[spin], s, mo_coeff[spin], optimize=True)

                #get maximum overlap
                idx = np.argsort(np.einsum('ij,ij->j', mom_s, mom_s))[::-1]
                mo_occ[spin][idx[:occupation_number]] = 1

            return mo_occ

        self.mom = get_orbital_occupation


    def execute(self, set=None):
        #main computation loop

        def compute_mo(f):
            #orthogonalise, solve and back-transform fock amtrix

            #orthogonalise Fock f->f' and solve
            fp = np.einsum('rp,xrs,sq->xpq', x, f, x, optimize=True )
            ep , cp = np.linalg.eigh(fp)

            #get the HOMO and LUMO energies
            homo = (ep[0][self.mol.nele[0]-1], ep[1][self.mol.nele[1]-1])
            lumo = (ep[0][self.mol.nele[0]],   ep[1][self.mol.nele[1]])

            #transform to ao basis
            c = np.einsum('pr,xrq->xpq', x, cp, optimize=True)

            return ep, c

        def get_energy(d, h_core, eri, ao, weights, cycle=0):
            #compute the 1e, coulomb and exchange-correlation energies

            one_electron_e = np.einsum('pq,pq->', d[0], h_core, optimize=True) + np.einsum('pq,pq->', d[1], h_core, optimize=True)
            j = np.einsum('rs,pqrs->pq', d[0]+d[1], eri, optimize=True)
            coulomb_e = np.einsum('pq,pq->', d[0]+d[1], j, optimize=True) * 0.5

            #evalute density over mesh
            rho = evaluate_rho_lda(d, ao, weights)*0.5
            
            #numerical electron count
            electrons = (np.sum(rho[0]*weights),np.sum(rho[1]*weights))
    
            #evaluate functional over mesh
            exc, vrho = lda_functional(self.x, self.c, rho)

            exc_e = np.sum(evaluate_exc(exc, rho, weights) )

            eSCF = one_electron_e + coulomb_e + exc_e

            out(self.mol.silent, [cycle, one_electron_e, coulomb_e, exc_e, electrons, eSCF],'scf')

            #evaluate potential
            vxc = evaluate_vxc(vrho, ao, weights)

            return vxc, eSCF

        #use a reduced version of Harpy's cython integrals
        from uks_aello import aello
        s, t, v, eri = aello(self.mol.atom, self.mol.orbital)

        #orthogonal transformation matrix
        from scipy.linalg import fractional_matrix_power as fractPow
        x = fractPow(s, -0.5).real

        #initial fock is core hamiltonian
        h_core = t + v

        #get alpha and beta electron counts
        paired_electrons = (np.sum(self.mol.nele) - self.mol.spin)//2
        self.mol.nele = [paired_electrons+self.mol.spin, paired_electrons]

        #get grid and weights
        grid, weights = GRID(self.mol.atom, self.mesh)

        #evaluate basis over grid
        ao = evaluate_atomic_orbital(self.mol.orbital, grid)

        #initial Fock guess
        f = (h_core, h_core)

        mo_energy, mo_coeff = compute_mo(f)

        #get occupied coefficients
        mo_occupation = self.mom(mo_energy, mo_coeff, s)

        #alpha and beta density matrix
        d = set if type(set) == np.ndarray else self.get_density_matrix(mo_coeff, mo_occupation)

        #unrestricted handler for closed closed shell
        mm = mo_coeff.copy()
        if (not self.open) and self.closed_shell_behavior != 'r' : 
            mo_coeff = self.closed_shell(mo_coeff)

        vxc, eSCF = get_energy(d, h_core, eri, ao, weights)

        last_cycle_energy = 0.0

        #diis initialisation
        if self.DIIS: diis = dii_subspace(self.diis_size)

        for cycle in range(1, self.cycles):
        	
            #construct Fock matrices 
            f = h_core + vxc +  np.einsum('rs,pqrs->pq', d[0]+d[1], eri, optimize=True)

            if (cycle != 0) and self.DIIS:
                f = diis.build(f, d, s, x)

            mo_energy, mo_coeff = compute_mo(f)

            #get occupied coefficients
            mo_occupation = self.mom(mo_energy, mo_coeff, s)

            d = self.get_density_matrix(mo_coeff, mo_occupation)

            vxc, eSCF = get_energy(d, h_core, eri, ao, weights, cycle)

            vector_norm = diis.norm if self.DIIS else ''
            out(self.mol.silent, [cycle, abs(eSCF - last_cycle_energy), vector_norm, self.get_spin_statistics(mo_coeff, mo_occupation, s) ], 'convergence')

            if abs(eSCF - last_cycle_energy) < self.tol: break

            last_cycle_energy = eSCF

        #final energies
        out(self.mol.silent, [eSCF, self.mol.nuclear_repulsion() ], 'final')

        #post SCF - mulliken and dipole and basis analysis
        out(self.mol.silent, [mo_energy, mo_coeff, self.mol.nele, d, s, self.mol], 'post')
        out(self.mol.silent, [self.mol], 'orbitals')


        #load cache with computed values
        self._cache['s'] = s ; self._cache['v'] = v ; self._cache['t'] = t ; self._cache['i'] = eri
        self._cache['f'] = f ; self._cache['d'] = d ; self._cache['c'] = mo_coeff ; self._cache['e'] = mo_energy ; self._cache['o'] = mo_occupation
        self._cache['g'] = grid ; self._cache['w'] = weights

        total_energy = eSCF + self.mol.nuclear_repulsion()
        print('Final Total Energy = {:>15.8f} Hartree'.format(total_energy))

        return total_energy
from __future__ import division
import numpy as np

class dii_subspace(object):
    #direct inversion of iterative subspace class

    def __init__(self, size):
        self.size = size
        self.fock_vector  = []
        self.error_vector = []
        self.nbf  = -1
        self.norm = 0.0

    def append(self, f, d, s, x):
        #update the subspaces respecting capacity of buffer

        if f.ndim == 3:
            #spin-polarised computation
            s = (s,s)
            g = []
            for ni in range(f.ndim-1):
                fds = np.einsum('im,mn,nj->ij',f[ni], d[ni], s[ni], optimize=True)
                g.append(np.einsum('mi,mn,nj->ij', x, (fds - fds.T), x, optimize=True))

            self.error_vector.append(np.vstack(g))
            self.fock_vector.append(f)
        else:
            #spin-unpolarised computation
            self.fock_vector.append(f) 
            fds = np.einsum('im,mn,nj->ij',f, d, s, optimize=True)
            self.error_vector.append(np.einsum('mi,mn,nj->ij', x, (fds - fds.T), x, optimize=True))

        #set nbf if not set
        if self.nbf == -1: self.nbf = f.shape[0]
        self.norm = np.linalg.norm(self.error_vector[-1])

        #check capacity
        if len(self.fock_vector) > self.size:
            del self.fock_vector[0]
            del self.error_vector[0]


    def build(self, f, d, s, x):
        #compute extrapolated Fock

        #update buffers
        self.append(f, d, s, x)

        #construct B matrix
        nSubSpace = len(self.fock_vector)
        b = -np.ones((nSubSpace+1,nSubSpace+1))
        b[:-1,:-1] = 0.0 ; b[-1,-1] = 0.0
        for i in range(nSubSpace):
            for j in range(nSubSpace):
                b[i,j] = b[j,i] = np.einsum('ij,ij->',self.error_vector[i], self.error_vector[j], optimize=True)


        #solve for weights
        residuals = np.zeros(nSubSpace+1)
        residuals[-1] = -1

        try:
            weights = np.linalg.solve(b, residuals)
        except np.linalg.LinAlgError as e:
            if 'Singular matrix' in str(e): exit('diis failed with singular matrix')

        #weights should sum to +1
        sum = np.sum(weights[:-1])
        assert np.isclose(sum, 1.0)

        #construct extrapolated Fock
        f = np.zeros_like(f, dtype='float')
        for i in range(nSubSpace):
            f += self.fock_vector[i] * weights[i]

        if f.ndim == 3:
            f = f.reshape(d.shape)

        return f

def out(silent, data, key):
    #print data to console

    symbol = lambda i: ['s','px','py','pz'][[[0,0,0],[1,0,0],[0,1,0],[0,0,1]].index(i)]

    if not silent and key == 'initial':

        units = {'bohr->angstrom' : 0.52917721092}

        mol = data[0]
        print()
        print('******************\n*   uks output   *\n******************')

        print('charge                  {:<2d}'.format(mol.charge))
        print('spin                    {:<2d}'.format(mol.spin))
        print('units                  ',mol.units)
        print('open shell             ',data[5], '     :  multiplicity ', int(2*(mol.spin/2)+1))

        print('\nbasis is               ', mol.name)
        print('analytic integration    aello cython - McMurchie-Davidson scheme')
        print('numerical integration   (Mura-Knowles, Lebedev)')
        print('                        radial prune is Aldrich-Treutler')
        print('                        Becke Partition scheme is Stratmann ')
        print('                        Radial adjust is Treutler')
        print('                        order: period 1 (10,11) and period 2 (15,15)')
        print('mesh                   ', data[4])
        print('functional             ', data[3].replace('&',','))
        print('\ndiis                   ', data[1],' : buffer size is ', data[2])

    if not silent and key == 'cycle':
        print('scf control             maximum cycles are ', data[0], '   :  convergence tolerance ', data[1])

    if not silent and key == 'scf':
        cycle, e1e, ej, ex, ne, et = data
        if cycle == 0:
            print('\n cycle     1 electron         coulomb         exchange          electrons           total')
            print('                                                                    S\u00B2      multiplicity            \u0394E         diis norm')
            print('----------------------------------------------------------------------------------------------------------------------------')
        print('   {:>2d}    {:>13.8f}     {:>12.8f}    {:>12.8f}    ({:>8.4f},{:>8.4f})  {:>12.8f} '.format(cycle, e1e, ej, ex, ne[0], ne[1], et))

    if not silent and key == 'convergence':
        cycle, delta, norm, spin_statistics = data
        ss, multiplicity = spin_statistics
        if cycle == 0:
            print()
        else:
            if type(norm) != str:
                print(' '*60, '     {:>6.3f}       {:>6.3f}         {:>12.6f}  {:>12.6f} '.format(ss, multiplicity, delta, norm))
            else: 
                print(' '*91, ' {:>12.6f}            '.format(delta))

    if not silent and key == 'final':
        print('\nnuclear repulsion   {:>15.10f}'.format(data[1]))
        print('total electronic    {:>15.10f}'.format(data[0]))
        print('\nfinal total energy  {:>15.10f}'.format(data[0]+data[1]))

    if not silent and key == 'post':
        #molecular orbitals for alpha and beta spins

        print('\nmolecular orbitals\n------------------')
        ep = data[0] ; cp = data[1] ; ne = data[2] ; d = data[3] ; s = data[4] ; mol = data[5] 
        spin = ['\u03B1', '\u03B2']
        for sp in range(2):
            print(spin[sp])
            print(np.array2string(ep[sp], floatmode='maxprec_equal', max_line_width=80, precision=5))
            print('HOMO {:<9.5f}  LUMO {:<9.5f} \n'.format(ep[sp][ne[sp]-1], ep[sp][ne[sp]]))

        #Lowdin population analysis
        print('\nLowdin populations\n--------------------')
        from scipy.linalg import fractional_matrix_power as fractPow
        u = fractPow(s, 0.5).real

        charge = np.zeros((2,mol.natm))
        for sp in range(2):
            pop = np.einsum('ij,jk,ki->i',u, d[sp], u, optimize=True)

            for i, a in enumerate(mol.atom):
                for j, o in enumerate(mol.orbital):
                    if o.atom.id == i: charge[sp][i] += pop[j]

            print(spin[sp])
            print(np.array2string(charge[sp], floatmode='maxprec_equal', max_line_width=80, precision=4))

        print(spin[0]+'+'+spin[1])
        resultant = charge[0] + charge[1]
        print(np.array2string(resultant, floatmode='maxprec_equal', max_line_width=80, precision=4))

        print('charge')
        for i,a in enumerate(mol.atom):
            resultant[i] = a.number - resultant[i]
        print(np.array2string(resultant, floatmode='maxprec_equal', max_line_width=80, precision=4), end='')
        print('        net = {:<5.2f}\n'.format(np.sum(resultant)))

        #alpha and beta dipole components
        from uks_aello import aello
        debyes = 2.541580253
        print('dipole momemts (Debye)\n----------------------')
        dipole = [np.array(aello(mol.atom, mol.orbital, 'dipole', d[0], [0,0,0])) * debyes * 0.5]
        dipole.append(np.array(aello(mol.atom, mol.orbital, 'dipole', d[1], [0,0,0])) * debyes * 0.5)
        print(' {:<3s}    x={:>8.4f} y={:>8.4f} z={:>8.4f}'.format(spin[0], dipole[0][0], dipole[0][1], dipole[0][2]))
        print(' {:<3s}    x={:>8.4f} y={:>8.4f} z={:>8.4f}'.format(spin[1], dipole[1][0], dipole[1][1], dipole[1][2]))
        resultant = dipole[0] + dipole[1]
        print(' {:<3s}    x={:>8.4f} y={:>8.4f} z={:>8.4f}'.format(spin[0]+'+'+spin[1], resultant[0], resultant[1], resultant[2]))

        print('\nGeometry\n-----------')
        for a in mol.atom:
            print('{:>2s}   {:>10.4f} {:>10.4f} {:>10.4f} '.format(a.symbol, a.center[0], a.center[1], a.center[2]))

    if key == 'orbitals':

        mol = data[0]
        if not silent:
            print('\nbasis analysis\n--------------')
            print('shell   l    pGTO    cGTO\n-------------------------')

        shells = 0; contracted = 0; primatives = 0
        basis_analysis = []

        for b in mol.orbital:
            basis_analysis.append([b.atom.id, np.sum(b.momenta),len(b.exponents), b.exponents[0]])

        nl = lambda x: (x+1)*(x+2)//2

        symbol = {0:'s',1:'p',2:'d'}
        from uks_mol import subshell_momenta

        mol_id = 0
        #group by atom center
        for a in mol.atom:
            atom_basis = [i for i in basis_analysis if i[0] == a.id]
            if not silent: print('{:<2s}   {:<2d}'.format(a.symbol, a.id))

            #group by momentum
            shell = [1, 2, 3]
            for l in range(3):
                momentum_basis = [i for i in atom_basis if i[1] == l]

                #polarized have multiple centers
                expected_contraction = len(momentum_basis)//nl(l)
                shell_basis = list(set(tuple(x) for x in momentum_basis))
                actual_contraction = len(shell_basis)

                contracted += len(shell_basis) * nl(l)

                #adjust contacted count for multiple center 
                delta_centers = expected_contraction - actual_contraction
                contracted += delta_centers * nl(l)

                for s in shell_basis:
                    centers = int([i[3] for i in momentum_basis].count(s[3])/nl(l))
                    if not silent: print('  {:<2d}    {:<2d}     {:<2d}      {:<1d}'.format(shells,s[1],s[2], centers))
                    shells += 1
                    primatives += s[2]*nl(l)

                    #construct the [n][l][m] designation of orbital
                    momentum = subshell_momenta[symbol[l]]
                    for center in range(centers):
                        for n in range(nl(l)):
                            i,j,k = momentum[n]
                            mol.orbital[mol_id].symbol = str(shell[l])+symbol[l]+'x'*i+'y'*j+'z'*k
                            mol_id += 1

                        shell[l] += 1
        if not silent:
            print('\n-------------------------')
            print('number of shells {:>12d}\nnumber of primative GTO {:>5d}\nnumber of contracted GTO {:>4d}'.format(shells, primatives, contracted))


def compute(ks, attribute='ionization potential',  dict=None):
    #perform some post scf computations

    from uks_uks import pKS

    eV = 27.211324570273

    if attribute == 'ionization potential':

        ground_state_energy = dict['eSCF']

        #get core-hole properties
        core_hole_coeff = ks.get('c')
        core_hole_occ   = ks.get('o')

        #core-hole at 1s alpha
        core_hole_occ[0][0] = 0

        #define new computation
        ch_ks = pKS(ks.mol, xc= ks.x + '&' + ks.c, mesh=ks.mesh)

        ch_ks.maximum_overlap_method(core_hole_coeff, core_hole_occ)
        core_hole_energy = ch_ks.execute()

        return (abs(ground_state_energy) - abs(core_hole_energy)) * eV

    elif attribute == 'excitation':
        import copy
        results = {} ; results['energy'] = []

        #get core-hole properties
        core_hole_coeff = ks.get('c')
        core_hole_occ   = ks.get('o')

        for level in dict['levels']:

            ground_state_energy = dict['eSCF']

            ch_occ   = core_hole_occ.copy()
            ch_coeff = core_hole_coeff.copy()

            ch_occ[1][level[0]] = 0
            ch_occ[1][level[1]] = 1

            d = ks.get_density_matrix(ch_coeff ,ch_occ)
            ks_ch = pKS(ks.mol, xc= ks.x + '&' + ks.c, mesh=ks.mesh)

            ks_ch.maximum_overlap_method(ch_coeff , ch_occ)
            e_excited = ks_ch.execute(set=d)

            results['energy'].append(((e_excited - ground_state_energy)*eV))

        return results['energy']


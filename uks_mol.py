from __future__ import division
import numpy as np

periodic_table = ['','H','He','Li','Be','B','C','N','O','F','Ne','Na']
subshell_momenta = {'s' : [(0,0,0)],
                    'p' : [(1,0,0),(0,1,0),(0,0,1)],
                    'd' : [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]}

class molecule(object):
    #molecule class - creates a molecule object
    '''
    properties - atom{symbol,number,center}
                 units   [string]
                 charge  [integer]
                 spin    [integer]
                 basis   [string] 
                 orbital{atom,primatives,momenta,exponents,coefficients,normals}

                 natm, norb, nele  
    methods    - nuclear_repulsion()               
    '''    
    def __init__(self, definition, units='angstrom', charge=0, spin=0, gto='sto-3g', path='', silent=False):
        
        class basis(object):

            def __init__(self, atom, primatives, momentum, exponents, coefficients, normals):

                self.atom = atom
                self.symbol = ''
                self.primatives = primatives
                self.momenta = np.array(momentum)
                self.exponents = np.array(exponents)
                self.coefficients = np.array(coefficients)
                self.normals = np.array(normals)

        class atoms(object):

            def __init__(self, id, symbol, center):

                self.id     = id
                self.symbol = symbol
                self.center = np.array(center, dtype=float)
                self.number = periodic_table.index(symbol)

        self.atom = []
        for n, a in enumerate(definition):
            self.atom.append(atoms(n, a[0],a[1]))

        self.units = units
        if units == 'angstrom':
            for i in self.atom:
                i.center /=  0.52917721092

        self.charge = charge
        self.spin = spin
        self.basis = basis
        self.name = gto

        self.orbital = []
        self.get_atom_basis(path)

        self.natm = len(self.atom)
        self.norb = len(self.orbital)
        self.nele = sum([i.number for i in self.atom]) - charge

        self.silent = silent

    def nuclear_repulsion(self):
            #nuclear repulsion 

        eNuc = 0.0
        natm = self.natm

        for i in range(natm):
            for j in range(i+1, natm):
                r = np.linalg.norm(self.atom[i].center - self.atom[j].center)
                eNuc += self.atom[i].number*self.atom[j].number / r

        return eNuc

    def get_atom_basis(self, path):
        #read atom blocks from BSE file

        def process_atom_block(block):
            #interpret BSE file data for atom (psi4 format)

            gto = {'m':[], 'p':[],'e':[],'c':[]}

            n = 0
            while True:

                m, p, _ = block[n].split()
                gto['m'].append(m.lower()) ; gto['p'].append(int(p))

                for i in range(gto['p'][-1]):
                    contractions = block[n+i+1].split()
                    gto['e'].append(float(contractions[0].replace('D','e')))
                    gto['c'].append(float(contractions[1].replace('D','e')))
                    if gto['m'][-1] == 'sp':
                        gto['c'].append(float(contractions[2].replace('D','e')))

                n += gto['p'][-1] + 1

                if n >= len(block)-1: return gto
                
        def create_basis_table(type_, gto, atom_table):
            #create instances of basis class

            ke = 0 ; kc = 0         
            for n,m in enumerate(gto['m']):

                if m != 'sp':
                    for L in subshell_momenta[m]:
                        atom_table.append(self.basis(type_, gto['p'][n], L, gto['e'][ke:ke+gto['p'][n]], gto['c'][ke:ke+gto['p'][n]], []))
                    ke += gto['p'][n] ; kc = ke

                else:
                    L = subshell_momenta['s'][0]
                    atom_table.append(self.basis(type_, gto['p'][n], L, gto['e'][ke:ke+gto['p'][n]], gto['c'][kc:kc+2*gto['p'][n]:2], []))
                    for L in subshell_momenta['p']:
                        atom_table.append(self.basis(type_, gto['p'][n], L, gto['e'][ke:ke+gto['p'][n]], gto['c'][kc+1:kc+2*gto['p'][n]:2], []))
                    ke += gto['p'][n] ; kc += 2*gto['p'][n]

            #sort orbitals by angular momentum
            atom_table.sort(key=lambda x: np.sum(x.momenta))

            return atom_table

        def normalisation_factor(orb):
            #compute the normalisation for contracted Gaussian

            #double factorial
            from scipy.special import factorial2 as df

            #princpal quantum number (n-1)
            n = np.sum(orb.momenta)

            #double factorial terms
            norm = np.zeros((orb.primatives))
 
            double_factorial_term = df(2*orb.momenta[0]-1) * df(2*orb.momenta[1]-1) * df(2*orb.momenta[2]-1)
            for i,p in enumerate(orb.exponents):
                norm[i] = pow(2*p/np.pi,0.75)*pow(4*p,n/2) / np.sqrt(double_factorial_term)

            #coefficients normalisation
            pf = pow(np.pi,1.5) * double_factorial_term / pow(2,n)

            s = 0.0
            for i, p in enumerate(orb.exponents):
                for j, q in enumerate(orb.exponents):
                    s += norm[i] * norm[j] * orb.coefficients[i] * orb.coefficients[j] / pow(p + q , n + 1.5)

            s *= pf
            s = 1/np.sqrt(s)
  
            orb.coefficients  *= s

            return orb.coefficients, norm


        def create_basis(atom_table, symbol, n):
            #assign atom basis to molecular atoms

            for i in atom_table:

                if i.atom == symbol:
                    i.coefficients, normal = normalisation_factor(i)
                    self.orbital.append(self.basis(self.atom[n], i.primatives, i.momenta, i.exponents, i.coefficients, normal)) 

        #get atomic species in order of increasing atomic number
        species = list(set([i.symbol for i in self.atom])) ; order = [periodic_table.index(i) for i in species]
        order = sorted(order)
        species = [periodic_table[i] for i in order]

        atom_table = []

        #read down file processing element list until empty
        fp = open(path, 'r')
        while fp:

            line = fp.readline().strip()

            if line.lower() in ['cartesian','spherical']: type_ = line

            #get the file block for atom type
            if line[:2].strip() == species[0]:
                atom_basis = []
                while line != '****':
                    line = fp.readline().strip()
                    atom_basis.append(line)

                #decode the atom block
                gto = process_atom_block(atom_basis)

                #produce basis for atom types
                atom_table = create_basis_table(species[0], gto, atom_table)

                species.pop(0)
                if species == []: break

        fp.close()

        #construct basis
        for n,i in enumerate(self.atom):
            create_basis(atom_table, i.symbol, n)

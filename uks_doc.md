1. **Defining The Molecule**

      mol = molecule([['O', (0.0, 0.0, 0.0)], 
                      ['H', (0,-0.757 ,0.587)], 
                      ['H', (0, 0.757 ,0.587)]], 
                      spin=1,
                      units='angstrom',
                      charge=-1,
                      gto='6-31g',
                      path='../harpy/basis/6-31g.gbf')

   The atom definitions should be item 1. They are a list of lists with each sub list being a string and a tuple. The string is an atomic element symbol (H-Na) and the    tuple is the x,y and z coordinates of the atom center.\
\
   The number of unpaired electrons is given by the *spin* keyword, the multiplicity is calculated as *spin*+1. For the *spin* to be odd the number of electrons must    also be odd, for the *spin* to be even the number of electrons must also be even. The program checks for this condition to be met.\
\
   The units are given by the *units* keyword and can be either 'angstrom' or 'bohr'. Internally coordinates are in 'bohr'\
\
   The overall charge on the molecule is given by the *charge* keyword as an integer.\
\
   The name of the basis is given by the *basis* keyword. This for display only as the actual basis used is determined by...\
\
   The Basis Set Exchange exchange file containing the Gaussian definitions for the basis is located by the *path* keyword. This should be a valid pathname.\
\
   The output (except final energies) can be suppressed by using the *silent* keyword. True will suppress output and (the default) False allows output.

2. **Molecular Properties.**\
   The atoms can be retrieved as, for example

      H1 = mol.atom[1]
        type = H1.symbol ; Z = H1.number ; x,y,z = H1.center
   Other supplied properties can be retrieved as, for example

        charge = mol.charge
   but note *gto* is returned as mol.name.

   The orbitals are found by

        Opx = mol.orbitals[2]
   with properties
   + Opx.symbol - orbital designation eg 1s, 2py
   + Opx.primatives   - the number of Gaussian primatives
   + Opx.coefficients - Gaussian primative coefficients (np array)
   + Opx.exponents - Gaussian prmative exponents (np array)
   + Opx.momenta - the components of the angular momentum of the orbital (np array)
   + Opx.normals - normalization factors for Gaussian primatives (np array)
   + Opx.atom - atom object for atom center of orbital

   Derived properties are
   + mol.natm - the number of atoms
   + mol.norb - the length of the basis
   + mol.nele - the number of electron
   + mol.nuclear_repulsion() - the nuclear repulsion energy

3. **Spin Polarized Parameters**\
   These are the parameters that define SCF procedure through the pKS class.

      ks = pKS(mol, 
              xc='LDA&VWN_RPA',
              cycles=50,
              tol=1e-8,
              mesh='close')

   The functional is given by the *xc* keyword. This is of the form '{exchange}&{correlation}', either {exchange} or {correlation} can be omitted. Currently 'LDA&VWN_RPA' is implemented.\
\
   The maximum number of iterations is given by the *cycles* keyword. Default is 50.\
\
The convergence tolerance is given by the *tol* keyword. Default is 1e-6.\
\
The integration grid spacing is given by the *mesh* keyword. This can be either 'coarse' or 'close' (see readme for details). Default is 'close'.\
\
Direct Inversion of Iterative Subspace is determined by the *diis* keyword which can be either True or False. Defaut is True.\
\
Size of diis buffer is given by the *diis_size* keyword. Default is 6.

4. **Spin Polarized Properties**\
   The *pKS* class has the following properties
   + ks.execute(set=) - this will initiate a spin polarized SCF computation. Returns SCF energy. The set keyword allows a density matrix pair to be supplied which will initialise the computation.
   + ks.open - boolean True for spin polarized, False for spin unpolarized. Returned not set.
   + ks.closed_shell_behavior - in event of a closed shell computation select either 'r' or 'u' to converge to a restricted or unrestricted solution. Position before calling kernel.
   + ks.get(key) - returns a computed quantity. *key* can be 's'-overlap, 't'-kinetic, 'v'-coulomb, 'i'-2 electron repulsion, 'f'-fock, 'd'-density, 'c'-MO coefficients, 'e'- MO energies, g'-grid coordinates and 'w'-grid weights. Note 'f', 'f', 'c' and 'e' return alpha and beta values.   We can find the difference between the exchange-correlation energy and the Hartree-Fock exchange energy as an example of ks.get() - see below.
   + ks.maximum_overlap_method(mo_coeff, mo_occ) - Imposes mo_occ occuaption numbers on the pending compuation. Position before calling kernel. See below for examples.

## Examples.
### The Difference between the exchange-correlation energy and the Hartree-Fock exchange energy.

    ks_energy = ks.execute()

    import numpy as np

    h_core = ks.get('t') + ks.get('v')                                      #core Hamiltonian

    eri = ks.get('i')  ;  d = ks.get('d')                                   #eri and density

    f =  np.einsum('pq,pq->', d[0]+d[1], h_core, optimize=True)             #one electron energy
    j =  np.einsum('rs,pqrs->pq', d[0]+d[1], eri, optimize=True) * 0.5      #coulomb
    f += np.einsum('pq,pq->', d[0]+d[1], j, optimize=True)                  #coulomb energy
    k = -np.einsum('rs,prqs->pq', d[0]+d[1], eri, optimize=True) * 0.25     #exchange
    f += np.einsum('pq,pq->', d[0]+d[1], k, optimize=True)                  #exchange energy
    f +=  mol.nuclear_repulsion()                                           #nuclear repulsion
    print(ks_energy - f)                                                    #delta
    >>>-0.03328974524120554

### Calculating the Ionization Potential using Maximum Overlap Method.

    from uks_util import compute

    mol = molecule([['C', (0,0,0)],
                    ['O', (0,0,1.147)]],
                    spin = 0,
                    units='angstrom',
                    charge=0,
                    gto='6-31g',
                    path='../harpy/basis/6-31g.gbf')

    ks = pKS(mol,
             xc='LDA&VWN_RPA')

    ks.closed_shell = 'r'
    ks_energy = ks.execute()

    mo_energy = ks.get('e')

    mol.silent=True
    ip =  compute(ks, 'ionization potential', dict={'eSCF': ks_energy} )

    print('Ionization Potential = {:<8.2f} eV'.format(ip))

    >>> Ionization Potential = 541.53   eV

Notes: *Path* is implemented so that if you have an existing application which already has basis sets then ther is no need to replicate those files. First we do a closed-shell calculation to get the ground state energy. The output is suppressed for the intermediate calculations by specifying 'mol.silent=True'. The ionization potential is then calculated by invoking the *compute* command from the uks_util module. The arguements for the call to *compute* (ks instance, method and a dictionary with a single entry 'eSCF' which is the ground state energy). The ionization potential is returned as a float.

### Excited states 
For C<sub>2</sub>H<sub>2</sub> in 6-31g basis the HOMO lies at (6+6+1+1)//2 ie orbital 7 (zero based) so we may calculate the k-edge of C1 as excitation 0->7

    from uks_util import compute
    import numpy as np

    def get_kedge_orbital(atom, coeff):
        #return the molecule orbital with s1 character for atom

        #get maximum ao value for each mo
        ao = np.argmax(np.absolute(coeff[0]), axis=0)

        for n,i in enumerate(ao):
            if (ks.mol.orbital[i].atom.id == atom) and (ks.mol.orbital[i].symbol == '1s'):
                return n

      mol = molecule([['C', (0,0,0)],
                      ['O', (0,0,1.147)]],
                      spin = 0,
                      units='angstrom',
                      charge=0,
                      gto='6-31g',
                      path='../harpy/basis/6-31g.gbf')

      mol.silent = True
      ks = pKS(mol, xc='LDA&VWN_RPA')
      ks.closed_shell = 'r'
      ks_energy = ks.execute()

      for atom in [0,1]:

          ao = get_kedge_orbital(atom, ks.get('c')) ; homo = ks.mol.nele[0]
          excitation = compute(ks, 'excitation', dict={'eSCF':ks_energy, 'levels':[(ao, homo)]})[0]
          print('K-Edge excitation for atom {:<1s} = {:<8.2f} eV'.format(mol.atom[atom].symbol, excitation))

    >>>Final Total Energy =   -112.63868065 Hartree
    >>>Final Total Energy =   -102.13586906 Hartree
    >>>K-Edge excitation for atom C = 285.80   eV
    >>>Final Total Energy =    -93.08753962 Hartree
    >>>K-Edge excitation for atom O = 532.01   eV

Notes: *compute* needs the ground state energy as before and a 'levels' keyword which is a list of tuples designating the excitations. *compute* returns a list of the excitation energies in eV. The routine *get_kedge_orbital(atom, ks)* finds the mo with 1s AO character for atom.

### Bond Orders

        import numpy as np
        geo = [[['C', (0, 0,  0.6146230)], ['C', (0, 0, -0.6146230)], ['H', (0, 0,  1.6926020)], ['H', (0, 0, -1.6926020)]],
               [['C', (0, 0, 0.6770110)], ['C', (0, 0, -0.6770110)], ['H', (0, 0.9347540, 1.2488340)],
                ['H', (0, -0.934754, 1.2488340)], ['H', (0, -0.934754, -1.2488340)], ['H', (0, 0.9347540, -1.2488340)]
               ]]



        mol = molecule( geo[0],
                        spin = 0,
                        units='angstrom',
                        charge=0,
                        gto='6-31g',
                        path='../harpy/basis/6-31g.gbf')

        mol.silent = False
        ks = pKS(mol, xc='LDA&VWN_RPA')
        ks.closed_shell = 'r'
        ks_energy = ks.execute()

        def mayer_bond_order(ks, bond):
            #get the Mayer Bond Order and valaence of atoms

            d = ks.get('d') ;  s = ks.get('s')
            p = np.einsum('xpr,rq->xpq', d, s, optimize=True)

            a = bond[0] ; b = bond[1]
           
            #get orbital set for each bond
            atom_orbitals = []
            for n, atom in enumerate(ks.mol.atom):
                atom_orbitals.append([k for k,i in enumerate(ks.mol.orbital) if i.atom.id == n])

            bos = []
            for n, atom in enumerate(ks.mol.atom):

                for i, bond in enumerate(ks.mol.atom):
                    if n == i: continue

                    s = 0
                    for mu in atom_orbitals[n]:
                        for nu in atom_orbitals[i]:
                            s += p[0][mu,nu]*p[0][nu,mu] + p[1][mu,nu]*p[1][nu,mu]  

                    bos.append(s)

            bond_order = 2*bos[int(a*ks.mol.natm + (b-1))]
            n = ks.mol.natm - 1
            valence = (2*np.sum(bos[a*n:a*n+n]), 2*np.sum(bos[b*n:b*n+n]))

            return bond_order, valence

    b, v = mayer_bond_order(ks, (0,1))

    print('bond C1-C2 has bond-order {:<4.2f}   C1 has valency {:<4.2f}   C2 has valency  {:<4.2f}  '.format(b, v[0], v[1]))

    >>>bond C1-C2 has bond-order 3.16   C1 has valency 4.07   C2 has valency  4.07  

and for C<sub>2</sub>H<sub>4</sub>

    >>>bond C1-C2 has bond-order 2.03   C1 has valency 3.89   C2 has valency  3.89  

Notes: The Mayer Bond order is described in this informative article [Bond order and valence indices: A personal account by I. Mayer](https://onlinelibrary.wiley.com/doi/10.1002/jcc.20494).
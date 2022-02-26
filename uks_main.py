from uks_mol import molecule
from uks_uks import pKS

if __name__ == '__main__':
    mol = molecule([['O', (0.0, 0.0, 0.0)], 
                    ['H', (0,-0.757 ,0.587)], 
                    ['H', (0, 0.757 ,0.587)]], 
                    spin=1,
                    units='angstrom',
                    charge=-1,
                    gto='6-31g',
                    path='6-31g.gbf')

    ks = pKS(mol, 
             xc='LDA&VWN_RPA',
             cycles=50,
             tol=1e-8,
             mesh='close')

    ks_energy = ks.execute()



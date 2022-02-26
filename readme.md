![image](https://user-images.githubusercontent.com/73105740/142720205-ddb5a6ad-4d6c-4a1a-8171-9ea1af5eece4.png)
### The kspy codes
  This is the fifth in the **kspy** series. The series is designed to elucidate various aspects of density functional theory. Each code is a self-contained example highlighting one particular feature, in this way the salient algorithms are not obscured with the code of a program structure. So far we have
  1.  **kspy_lda**\
     Introduced an atom-centered integration grid, direct-inversion of iterative subspace, molecular integrals, closed-shell self-consistent calculation and Local Spin Density Approximation. The LSDA functional illustrated is Slater exchange and VWN3 correlation.
  2.  **kspy_gga**\
     Adds routines to compute Gaussian gradients over an integration grid and implements a Generalized Gradient Approximation functional. The functionals implemented are PBE and B3LYP. A hybrid PBE0 functional is also implemented to illustrate hybrid-GGA. 
  3.  **kspy_mgga**\
     Adds routines to compute Gaussian second derivatives over an integration grid and implements a meta-Generalized Gradient Approximation (mGGA). The functional implemented is MS (made-simple).
  4.  **kspy-tddft**\
     This implements Time-Dependent DFT theory in the Tamm-Dancoff Approximation (TDA) by solving the Casida equation. Linear response transition properties such as dipoles, oscillator and rotary strengths are computed together with natural transition orbitals. The transition energy spectrum can be plotted with broadening. Real-time TDA is implemented via density matrix propogation using the Magnus second-order propogator. Absorption spectra can be plotted using a Pade approximants approach. The functionals used are Slater exchange and RPA correlation.
  5.  **kspy-uks**\
     This is finally a 'program'. It illustrates creating a basis from a Basis Set Exchange (BSE), modification to DIIS for alpha and beta spin matrices, the solution of the spin polarized Kohn-Sham equations. The functional used is the Slater exchange with VWN(3) RPA correlation. The [Maximum Overlap Method](https://rsc.anu.edu.au/~pgill/papers/118MOM.pdf) is also implemented.

### Density Function Theory Program - kspy-sp
  The program **kspy-lda** showed how to implement a Mura-Knowles/Lebedev grid and use a LDA functional. This program implements a Generalized Gradient Approximation (GGA) based functional - the VWN_RPA. The previous programs were in the spin unpolarized (restricted) formalism, whereas here we go to the full spin polarized (unrestricted) Kohn-Sham equations. 

### The Grid
  The grid is the same as **kspy-lda**. The radial grid is a Mura-Knowles radial grid [ME Mura and PJ Knowles 'Improved radial grids for quadrature in density-functional calculations' JCP 104, 9848 (1996); DOI:10.1063/1.471749](https://aip.scitation.org/doi/10.1063/1.471749). The 'coarse' angular grid is of Lebedev orders (11, 15) for period 1 and period 2 respectively. This translates into 50 and 86 points respectively arranged on a spherical shell (VI Lebedev, and DN Laikov, Doklady Mathematics, 'A Quadrature formula for the sphere of the 131st algebraic order of accuracy' Vol. 59, No. 3, (1999)). There are various sources for this data given in the external links of the wikipedia article on Lebedev integration.  A pruning scheme is employed to systematically reduce the number of angular points in regions where dense angular quadrature is not necessary, such as near the nuclei where the charge density is approximately spherically symmetric and at long distance from the nucleus. The pruning scheme I employed was the Treutler-Aldrich scheme [O Treutler and R Ahlrich, 'Efficient molecular numerical integration schemes' JCP 102, 346 (1995); DOI:10.1063/1.469408](https://aip.scitation.org/doi/pdf/10.1063/1.469408). The partitioning of the atomic centered grids to a molecular grid follows a Becke scheme after Stratmann [RE Stratmann, GE Scuseria and MJ Frisch, 'Achieving Linear scaling in exchange-correlation density functional quadratures' CPL 257, 3-4 (1996); DOI:10.1016/009-2614(96)00600-8](https://www.sciencedirect.com/science/article/abs/pii/0009261496006008?via%3Dihub). Finally I have implemented a final radius adjustment during the partition (Becke suggests doing this) using the Bragg radius. A second 'close' grid is also included which is a (50, 75) radial and (29, 29) angular, the latter representing 302 points on each shell. The grid routines are in ks_grid.py.

### The HF Integrals
  To get the DFT SCF started we need an initial density. To do this I use a HF overlap matrix S, and an initial Fock matrix composed of the sum of the 1-electron kinetic and coulomb integrals (core Hamiltonian - T+V). This Fock is then orthogonalised (F<sup>'</sup>) as (S<sup>-0.5</sup>)<sup>T</sup>FS<sup>-0.5</sup>, eigensolve the resulting orthogonal Fock for orbital coefficients C orthogonal, transform back to atomic basis as S<sup>-0.5</sup>C<sup>'</sup>, use resulting ao coefficients to compute a density matrix D<sub>&mu;&nu;</sub> = c<sub>&mu;i</sub>c<sub>i&nu;</sub> where i is over occupied orbitals. This initial density can be used with initial Fock and 2-electron repulsion integrals to form the coulomb integral J and the exchange integral K. To get these integrals I've used a modified version of Harpy's Cython integral package *aello*. I've removed the electric field integrals and additionally the 2-electron repulsions are returned as a tensor rather than a linear array. These are in uks_aello.pyx.

### Molecule and Basis Sets
  The molecule definition is contained in a *mol* object which is itself comprised of objects from an atom class and a basis class. Each instance of the atom class contains an id, atom symbol, atomic number and the coordinates of the atom center (array[3]). Each instance of the basis class contains the atom object the Gaussian is centered on, the momentum(array[3]), the exponents (array[primatives]), the coefficients (array[primatives]) and the normalisation (array[primatives]). The momenta are given as eg s [0,0,0] px [1,0,0] py [0,1,0] and pz [0,0,1]. The mol class reads a BSE (psi4 format) basis file and creates a mol class instance which includes an array of basis objects. The basis file used is specified by the *path* keyword so existing basis set file can be used, however I have included 6-31g in case you don't already have a basis set file.

### The Functional
  We are using the GGA (Generalized Gradient Approximation) with the VWN RPA functional - defined [here](https://escholarship.org/content/qt23j4q7zm/qt23j4q7zm.pdf?t=obc5l4). The functional has been implemented in a way that is compatible with **libxc**. That is to say the LDA energy is included and rather than the pure gradients a contracted version (&sigma;) are passed to the functional evaluator. &sigma; is defined as &Delta;&rho;<sub>i</sub>.&Delta;&rho;<sub>i</sub>, details are given in the **libxc** documentation (https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/). The functional code and gradient evaluation routines are in uks_xc.py. 

### The Density
  The density is determined by first evaluating each Gaussian orbital over the grid (&rho;<sub>g</sub>), then assembling those into a complete basis over the grid (&rho;<sub>b</sub> = &Sigma;&rho;<sub>g</sub>) and finally computing &rho; as the product  &rho;<sub>b</sub>*&rho;<sub>b</sub>D. 

### Convergence
  The diis scheme is a modified version for spin polarized systems of that in kspy_lda. This is implemented as a class in uks_util.

### Output
    ******************
    *   uks output   *
    ******************
    charge                  -1
    spin                    1 
    units                   angstrom
    open shell              True         multiplicity  2

    basis is                6-31g
    analytic integration    aello cython - McMurchie-Davidson scheme
    numerical integration   (Mura-Knowles, Lebedev)
                            radial prune is Aldrich-Treutler
                            Becke Partition scheme is Stratmann 
                            Radial adjust is Treutler
                            order: period 1 (10,11) and period 2 (15,15)
    mesh                    close
    functional              LDA,VWN_RPA

    diis                    True   buffer size is  6
    scf control             maximum cycles are  50         convergence tolerance  1e-08

     cycle     1 electron         coulomb         exchange          electrons           total
                                                                    S²      multiplicity            ΔE         diis norm
    ----------------------------------------------------------------------------------------------------------------------------
    0     -139.68763205     73.83276889    -12.12070629    (  6.0000,  5.0000)  -77.97556946 
    1     -104.88495728     33.78263452     -7.50592719    (  6.0000,  5.0000)  -78.60824995 
                                                                   1.758        2.834            78.608250      3.492575 
    2     -124.14047494     48.17056327     -8.98224474    (  6.0000,  5.0000)  -84.95215641 
                                                                   0.756        2.006             6.343906      3.969463 
    3     -123.95728101     47.97673872     -8.98375811    (  6.0000,  5.0000)  -84.96430039 
                                                                   0.752        2.002             0.012144      0.768064 
    4     -124.35930767     48.39288095     -9.02848055    (  6.0000,  5.0000)  -84.99490727 
                                                                   0.752        2.002             0.030607      0.708567 
    5     -126.01723618     50.19490458     -9.23473750    (  6.0000,  5.0000)  -85.05706910 
                                                                   0.751        2.001             0.062162      0.587773 
    6     -126.16333066     50.35883931     -9.25323844    (  6.0000,  5.0000)  -85.05772979 
                                                                   0.751        2.001             0.000661      0.061987 
    7     -126.20487852     50.40580709     -9.25868758    (  6.0000,  5.0000)  -85.05775901 
                                                                   0.751        2.001             0.000029      0.012652 
    8     -126.20273484     50.40338613     -9.25841039    (  6.0000,  5.0000)  -85.05775910 
                                                                   0.751        2.001             0.000000      0.000660 
    9     -126.20282114     50.40348329     -9.25842125    (  6.0000,  5.0000)  -85.05775910 
                                                                   0.751        2.001             0.000000      0.000032 

    nuclear repulsion      9.1882584177
    total electronic     -85.0577591019

    final total energy   -75.8695006841


### Test
  The python code for the functional used here was copied into a pyscf eval_xc subroutine (user defined functional) and run with the same grid. The pyscf value is -75.8695006841355 Hartree <sup>*</sup>, this agrees better than the convergence criteria of 1e-8. The value using libxc from pyscf is -75.8695006841354 Hartree. 
  
  **<sup>*</sup>** the basis definition with pyscf (nwchem) differs slightly from the basis here (psi4).

### Installation
  Copy all files to a directory, run

    python3 setup.py build_ext --inplace install --user

then 

    python3 uks_main.py

### Usage
   The file uks_main contains an example of how to run the code, further (more detailed) documentation can be found in uks_doc.md.
   
### Post SCF
    molecular orbitals
    ------------------
    α
    [-16.19913  -0.18130   0.19314   0.35737   0.42795   0.53869   0.62960
       1.39344   1.45493   1.46704   1.56804   1.61609   2.00151]
    HOMO 0.53869    LUMO 0.62960   

    β
    [-16.20010  -0.18675   0.18614   0.35213   0.42468   0.51638   0.60809
       1.38662   1.45042   1.46081   1.56079   1.60794   1.99503]
    HOMO 0.42468    LUMO 0.51638   


    Lowdin populations
    --------------------
    α
    [4.4638 0.7681 0.7681]
    β
    [4.1326 0.4337 0.4337]
    α+β
    [8.5964 1.2018 1.2018]
    charge
    [-0.5964 -0.2018 -0.2018]        net = -1.00

    dipole momemts (Debye)
    ----------------------
     α      x= -0.0000   y= -0.0000  z=  -2.3595
     β      x=  0.0000   y= -0.0000  z=   0.7741
     α+β    x=  0.0000   y= -0.0000  z=  -1.5854

Finally the geometry and basis are reported (pGTO are the number of primatives and cGTO are the number of centers)

    Geometry
    -----------
     O       0.0000     0.0000     0.0000 
     H       0.0000    -1.4305     1.1093 
     H       0.0000     1.4305     1.1093 

    basis analysis
    --------------
    shell   l    pGTO    cGTO
    -------------------------
    O    0 
      0     0      6       1
      1     0      3       1
      2     0      1       1
      3     1      1       1
      4     1      3       1
    H    1 
      5     0      3       1
      6     0      1       1
    H    2 
      7     0      3       1
      8     0      1       1

    -------------------------
    number of shells            9
    number of primative GTO    30
    number of contracted GTO   13

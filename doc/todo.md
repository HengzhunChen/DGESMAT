TODO List
=================================================================

General
-----------------------------------------------------------------
- Implement ionic dynamics

- Implement multiple spin component wave function, see also Spinor class

- Add support for isUseVLocal for HGH type pseudopotential

- Add support for more type of exchange-correlation functional, currently 
  not hybrid type XC functional (e.g. HSE06) is supported. Corresponding 
  VDW energy and force is also need to supported. See also LibXC folder.

- Test smearing option, currently have NOT suitable test example
 
- In CalculateNonlocalPP(), move spherical harmonic function into fodler 
  Utilities to simplify the code

- Add plot function in Utilities folder

- Fix the cutoff in ReadUPF() since they should be given by a table 
  according to the element type but not use the same value


PWDFT
-----------------------------------------------------------------
- Implement Hybrid


DGDFT
------------------------------------------------------------------
- Refactor the EigenSolverKS and move the eigen-decomposition solver to 
  utilities folder so that DG_Sovler can also use them

- Complete CheFSI for DG_Solver
 
- Implement PEXSI as option of DG_Solver

- Add support for user option isUseVLocal for DG method

- Add support for user option isCalculateForceEachSCF 

- Add support for user option isCalculateAPosterioriEachSCF

- Test the calculation of Efree in finite temperature, currently have NOT 
  suitable test example
 
- Fix the mtimes() of HamiltonianDG, currently have not utilized the 
  sparse structure of DG Hamiltonian matrix than just treat it as a sparse 
  matrix. NOTE: this refactor may not reduce run time comparing with eigs 
  for sparse matrix in matlab, need to check by some experiment.

- Find the convergence problem in calculating force together with the 
  ONCV pseudopotential

- Add more suitable examples to test dgdft in test folder


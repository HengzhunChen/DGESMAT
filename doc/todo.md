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
 
- In CalculateNonlocalPP(), move spherical harmonic function into folder 
  Utilities to simplify the code

- Add plot function in Utilities folder

- Fix the cutoff in ReadUPF() since they should be given by a table 
  according to the element type but not use the same value


DGDFT
------------------------------------------------------------------
- Add support for user option isUseVLocal for DG method

- Add support for user option isCalculateAPosterioriEachSCF

- Test the calculation of Efree in finite temperature, currently have NOT 
  suitable test example
 
- Add more suitable examples to test dgdft in test folder


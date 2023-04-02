TODO List
=================================================================

General
-----------------------------------------------------------------
- Implement multiple spin component wave function, see also Spinor class

- Add support for more type of exchange-correlation functional, currently 
  not hybrid type XC functional (e.g. HSE06) is supported. Corresponding 
  VDW energy and force is also need to supported. See also LibXC folder.

- Test smearing option and clean relative code. 
 
- In CalculateNonlocalPP(), move spherical harmonic function into folder 
  Utilities to simplify the code

- Add plot function in Utilities folder

- Consider data structure of potential in @PeriodTable, whether to keep 
  the derivative of potential together with idx and val. 

- Check the influence of the adjustment to pseudocharge in function
  CalculatePseudoCharge


DGDFT
------------------------------------------------------------------
- Add support for user option isCalculateAPosterioriEachSCF

- Test the calculation of Efree in finite temperature


Remark
===================================================================

Some features of DGESMAT differ with ScalES
-------------------------------------------------------------------
- Fourier interpolation from coarse grid to fine grid, see also 
  @Fourier, @Spinor, @HamiltonianKS

- Some default values are different, see also ESDFReadInput.m

- Cutoff and RGaussian may be different, see also ReadUPF.m

- Add many warning and error to provide guidance to user

- A check for multi-defined atom in extended element, see also
  setup_element.m

- Implement some auxiliary functions to make the DFT calculation and 
  post-processing easier to manipulate.

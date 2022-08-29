function define_global(outFile, debugFile)
% DEFINE_GLOBAL define global variables
%
%    define_global(outFile, debugFile) initializes some global variables 
%    and is used as a reference for the usage of the global variables.
%
%    See also dgdft_main, pwdft_main.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;
esdfParam = ESDFInputParam();
% used in
% ESDF/
%    ESDFReadInput()
%    ESDFPrintInput()
% EigenSolverKS/
%    PPCGSolve()
% HamiltonianDG/
%    CalculateAtomDensity()
%    Setup()
% HamiltonianKS/
%    CalculateAtomDensity()
%    CalculateForce()
%    CalculateIonSelfEnergyAndForce()
%    CalculatePseudoPotential()
%    CalculateVdwEnergyAndForce()
%    CalculateVtot()
%    Setup()
% SCF/
%    CalculateEnergy()
%    Iterate()
%    Setup()
% SCFDG/
%    InnerIterate()
%    Iterate()
%    Setup()
%    UpdateElemLocalPotential()
% PeriodTable/
%    SelfIonInteraction()
%    Setup()
% pwdft_main()
% dgdft_main()
% setup_element()


% output files
global outFid
outFid = fopen(outFile, 'w');
% debug output file
global debugFid
if debugFile ~= ""
    debugFid = fopen(debugFile, 'w');
else
    debugFid = [];
end
% used in
% pwdft_main()
% dgdft_main()
% InfoPrint()
% PrintBlock()


%
% global variables used in esdfParam to read input
%

% global variables deal with conversion of different physical units
global phy_var phy_name phy_unit;
phy_var = "";
phy_name = "";
phy_unit = [];
% used in esdf_unit(), esdf_convfac()

global kw_label lineList tokenList;
% kw_label is the list of all ESDF parameters
% lineList is the line list read from input file, used for extracting
% tokenList and read data for block type.
% tokenList is the list of tokens, each row is label and its value, not
% including the block type labels.
kw_label = "";
lineList = "";
tokenList = "";
% used in esdf_init(), esdf_key(), esdf_block(), esdf_get()


end
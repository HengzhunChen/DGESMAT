function HamKS = Setup(HamKS)
% HAMILTONIANKS/SETUP initializes HamiltonianKS object HamKS with data 
%    from global variable esdfParam.
% 
%    See also HamiltonianKS, Domain, Atom, ESDFInputParam.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;

% FIXME Hard coded
HamKS.numSpin = 2;

HamKS.numExtraState = esdfParam.basic.numExtraState;


ntotFine = HamKS.domain.NumGridTotalFine();
dim = dimDef();

HamKS.density = zeros(ntotFine, 1);
HamKS.gradDensity = cell(dim, 1);
for d = 1 : dim
    HamKS.gradDensity{d} = zeros(ntotFine, 1);
end


HamKS.pseudoCharge = zeros(ntotFine, 1);
if esdfParam.userOption.general.isUseVLocal
    HamKS.vLocalSR = zeros(ntotFine, 1);
    HamKS.gaussianCharge = zeros(ntotFine, 1);
end
HamKS.vext  = zeros(ntotFine, 1);
HamKS.vhart = zeros(ntotFine, 1);
HamKS.vtot  = zeros(ntotFine, 1);
HamKS.epsxc = zeros(ntotFine, 1);
HamKS.vxc   = zeros(ntotFine, 1);

% initialize the XC functionals
HamKS.isHybrid = false;
HamKS.XCType = esdfParam.basic.XCType;
if contains(HamKS.XCType, "HYB")
    HamKS.isHybrid = true;
end


% Set up wavefunction filter options, useful for CheFSI in PWDFT
if esdfParam.basic.PWSolver == "CheFSI"
    HamKS.ecutFilter.isApplyFilter = 1;
    HamKS.ecutFilter.wfnCutoff = esdfParam.basic.ecutWavefunction;
else
    HamKS.ecutFilter.isApplyFilter = 0;
    HamKS.ecutFilter.wfnCutoff = esdfParam.basic.ecutWavefunction;
end    


% parameters for hybrid
HamKS.hybrid.DFType              = esdfParam.hybrid.DFType;
HamKS.hybrid.DFKmeansWFType      = esdfParam.hybrid.DFKmeansWFType;
HamKS.hybrid.DFKmeansWFAlpha     = esdfParam.hybrid.DFKmeansWFAlpha;
HamKS.hybrid.DFKmeansTolerance   = esdfParam.hybrid.DFKmeansTolerance;
HamKS.hybrid.DFKmeansMaxIter     = esdfParam.hybrid.DFKmeansMaxIter;
HamKS.hybrid.DFNumMu             = esdfParam.hybrid.DFNumMu;
HamKS.hybrid.DFNumGaussianRandom = esdfParam.hybrid.DFNumGaussianRandom;
HamKS.hybrid.DFTolerance         = esdfParam.hybrid.DFTolerance;
HamKS.hybrid.exxDivergenceType   = esdfParam.hybrid.exxDivergenceType;


end
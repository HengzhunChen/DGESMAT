function HamKS = Setup(HamKS, esdfParam)
% HAMILTONIANKS/SETUP initializes HamiltonianKS object HamKS with data 
%    from global variable esdfParam.
% 
%    See also HamiltonianKS, Domain, Atom, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


HamKS.isDGDFT = esdfParam.isDGDFT;

HamKS.ecutWavefunction = esdfParam.basic.ecutWavefunction;

% FIXME Hard coded
HamKS.numSpin = 2;

HamKS.numExtraState = esdfParam.basic.numExtraState;
HamKS.numExtraElectron = esdfParam.basic.extraElectron;

ntotFine = HamKS.domain.NumGridTotalFine();
dim = dimDef();

HamKS.density = zeros(ntotFine, 1);
HamKS.gradDensity = cell(dim, 1);
for d = 1 : dim
    HamKS.gradDensity{d} = zeros(ntotFine, 1);
end


HamKS.pseudoCharge = zeros(ntotFine, 1);
HamKS.vLocalSR     = zeros(ntotFine, 1);
HamKS.vext  = zeros(ntotFine, 1);
HamKS.vhart = zeros(ntotFine, 1);
HamKS.vtot  = zeros(ntotFine, 1);
HamKS.epsxc = zeros(ntotFine, 1);
HamKS.vxc   = zeros(ntotFine, 1);

HamKS.XCType     = esdfParam.basic.XCType;
HamKS.pseudoType = esdfParam.basic.pseudoType;
HamKS.VDWType    = esdfParam.basic.VDWType;


% Set up wavefunction filter options, useful for CheFSI in PWDFT
if esdfParam.basic.PWSolver == "CheFSI"
    HamKS.ecutFilter.isApplyFilter = 1;
    HamKS.ecutFilter.wfnCutoff = esdfParam.basic.ecutWavefunction;
else
    HamKS.ecutFilter.isApplyFilter = 0;
    HamKS.ecutFilter.wfnCutoff = esdfParam.basic.ecutWavefunction;
end    


end
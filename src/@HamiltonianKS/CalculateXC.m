function [Exc, epsxc, vxc] = CalculateXC(HamKS)
% HAMILTONIANKS/CALCULATEXC calculates exchange correlation energy and
%    potential.
%
%    See also HamiltonianKS, xcRef.    

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


ntotFine = HamKS.domain.NumGridTotalFine();
vol = HamKS.domain.Volume();

F = HamKS.fft;

% cutoff
epsRho = 1e-8;
epsGRho = 1e-10;

XCType = HamKS.XCType;

if contains(XCType, "LDA")
    [epsxc, vxc] = xc_lda_exc_vxc(XCType, HamKS.density);
    
    % modify "bad points"
    idxzero = HamKS.density < epsRho;
    epsxc(idxzero) = 0;
    vxc(idxzero) = 0;
    
% ---------------- end of XC_FAMILY_LDA -------------------------- 
    
elseif contains(XCType, "GGA")    
    gradDensityX = HamKS.gradDensity{1};
    gradDensityY = HamKS.gradDensity{2};
    gradDensityZ = HamKS.gradDensity{3};

    gradDensity2 = gradDensityX.^2 + gradDensityY.^2 + gradDensityZ.^2;

    density = HamKS.density;
    
    [epsxc, vxc1, vxc2] = xc_gga_exc_vxc(XCType, density, gradDensity2);

    % modify "bad points"
    idxz = density < epsRho | gradDensity2 < epsGRho;
    epsxc(idxz) = 0;
    vxc1(idxz)  = 0;
    vxc2(idxz)  = 0;
    
    vxc = vxc1;
    for d = 1 : dimDef()
        gradDensityd = HamKS.gradDensity{d};
        x = gradDensityd .* 2 .* vxc2;
        y = F * x;
        ik = F.ikFine{d};
        idxnz = F.gkkFine ~= 0;
        y(~idxnz) = 0;
        y(idxnz) = y(idxnz) .* ik(idxnz);
        x = F' * y;
        vxc = vxc - real(x);
    end
        
% -------------------- end of XC_FAMILY_GGA --------------------------
    
else
    error('Unsupported XC family!');
end


% ----------- compute total exchange-correlation energy ----------------

Exc = sum(HamKS.density .* epsxc) * vol / ntotFine;


end
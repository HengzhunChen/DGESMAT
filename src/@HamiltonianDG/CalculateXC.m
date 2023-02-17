function [Exc, epsxc, vxc] = CalculateXC(HamDG)
% HAMILTONIANDG/CALCULATEXC calculates exchange correlation energy and
%    potential.
%
%    [Exc, epsxc, vxc] = CalculateXC(HamDG) returns exchange correlation
%    energy Exc and potential over each element and stores as cells 
%    epsxc, vxc. 
%
%    See also HamiltonianDG, xcRef.    

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


F = HamDG.fft;

Exc = 0;
numElem = HamDG.numElem;
numElemTotal = prod(numElem);
epsxc = cell(numElemTotal, 1);
vxc = cell(numElemTotal, 1);

% Cutoff of the XC potential. Important for SCF convergence for GGA and
% above.
epsRho = 1e-8;
epsGRho = 1e-8;

XCType = HamDG.XCType;

if contains(XCType, "LDA")
    for elemIdx = 1 : numElemTotal
        localRho = HamDG.density{elemIdx};
        
        [localEpsxc, localVxc] = xc_lda_exc_vxc(XCType, localRho);

        % modify "bad points"
        idxzero = localRho < epsRho;
        localEpsxc(idxzero) = 0;
        localVxc(idxzero) = 0;
        
        epsxc{elemIdx} = localEpsxc;
        vxc{elemIdx} = localVxc;
        Exc = Exc + sum(localRho .* localEpsxc);
    end
% -------------------- end of XC_FAMILY_LDA -----------------------------
    
elseif contains(XCType, "GGA")
    vxc2 = cell(numElemTotal, 1);
    for elemIdx = 1 : numElemTotal
        localRho = HamDG.density{elemIdx};
        gradDensityX = HamDG.gradDensity{1}{elemIdx}; 
        gradDensityY = HamDG.gradDensity{2}{elemIdx}; 
        gradDensityZ = HamDG.gradDensity{3}{elemIdx};

        localGRho2 = gradDensityX.^2 + gradDensityY.^2 + gradDensityZ.^2;

        [localEpsxc, localVxc, localVxc2] = ...
            xc_gga_exc_vxc(XCType, localRho, localGRho2);

        % modify "bad points"
        idxz = localRho < epsRho | localGRho2 < epsGRho;
        localEpsxc(idxz) = 0;
        localVxc(idxz) = 0;
        localVxc2(idxz) = 0;

        epsxc{elemIdx} = localEpsxc;
        vxc{elemIdx} = localVxc;
        vxc2{elemIdx} = localVxc2;

        Exc = Exc + sum(localRho .* localEpsxc);
    end

    for d = 1 : dimDef()
        gradDensityVxc2 = cell(numElemTotal, 1);
        for elemIdx = 1 : numElemTotal
            localVxc2 = vxc2{elemIdx};
            gradDensityd = HamDG.gradDensity{d}{elemIdx};
            
            localGradDensityVxc2 = gradDensityd .* 2 .* localVxc2;
            gradDensityVxc2{elemIdx} = localGradDensityVxc2;
        end

        tempVec = ElemVecToGlobal(gradDensityVxc2, HamDG.domain.numGridFine, numElem);

        Y = F * tempVec;
        ik = F.ikFine{d};
        idxnz = F.gkkFine ~= 0;
        Y(idxnz) = Y(idxnz) .* ik(idxnz);
        Y(~idxnz) = 0;
        X = real(F' * Y);

        gradGradDensityVxc2 = GlobalVecToElem(X, HamDG.domain.numGridFine, numElem);

        for elemIdx = 1 : numElemTotal
            localVxc = vxc{elemIdx};
            localGradGradDenistyVxc2 = gradGradDensityVxc2{elemIdx};
            localVxc = localVxc - localGradGradDenistyVxc2;
            vxc{elemIdx} = localVxc;
        end
    end
% -------------------- end of XC_FAMILY_GGA ----------------------------

else
    error('Unsupported XC family!');
end

Exc = Exc * HamDG.domain.Volume() / HamDG.domain.NumGridTotalFine();

end
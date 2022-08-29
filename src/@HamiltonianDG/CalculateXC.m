function [Exc, epsxc, vxc] = CalculateXC(HamDG)
% HAMILTONIANDG/CALCULATEXC calculates exchange correlation energy and
%    potential.
%
%    [Exc, epsxc, vxc] = CalculateXC(HamDG) returns exchange correlation
%    energy Exc and potential over each element and stores as cells 
%    epsxc, vxc. 
%
%    See also HamiltonianDG, xcRef.    

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


F = HamDG.fft;

Exc = 0;
numElem = HamDG.numElem;
epsxc = cell(numElem);
vxc = cell(numElem);

% Cutoff of the XC potential. Important for SCF convergence for GGA and
% above.
epsRho = 1e-8;
epsGRho = 1e-8;

XCType = HamDG.XCType;

if contains(XCType, "LDA")
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                localRho = HamDG.density{i, j, k};
                
                [localEpsxc, localVxc] = xc_lda_exc_vxc(XCType, localRho);

                % modify "bad points"
                idxzero = localRho < epsRho;
                localEpsxc(idxzero) = 0;
                localVxc(idxzero) = 0;
                
                epsxc{i, j, k} = localEpsxc;
                vxc{i, j, k} = localVxc;
                Exc = Exc + sum(localRho .* localEpsxc);
            end
        end
    end
% -------------------- end of XC_FAMILY_LDA -----------------------------
    
elseif contains(XCType, "GGA")
    vxc2 = cell(numElem);
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                localRho = HamDG.density{i, j, k};
                gradDensityX = HamDG.gradDensity{1}{i, j, k}; 
                gradDensityY = HamDG.gradDensity{2}{i, j, k}; 
                gradDensityZ = HamDG.gradDensity{3}{i, j, k};

                localGRho2 = gradDensityX.^2 + gradDensityY.^2 + gradDensityZ.^2;

                [localEpsxc, localVxc, localVxc2] = ...
                    xc_gga_exc_vxc(XCType, localRho, localGRho2);

                % modify "bad points"
                idxz = localRho < epsRho | localGRho2 < epsGRho;
                localEpsxc(idxz) = 0;
                localVxc(idxz) = 0;
                localVxc2(idxz) = 0;

                epsxc{i, j, k} = localEpsxc;
                vxc{i, j, k} = localVxc;
                vxc2{i, j, k} = localVxc2;

                Exc = Exc + sum(localRho .* localEpsxc);
            end
        end
    end

    for d = 1 : dimDef()
        gradDensityVxc2 = cell(numElem);
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    localVxc2 = vxc2{i, j, k};
                    gradDensityd = HamDG.gradDensity{d}{i, j, k};
                    
                    localGradDensityVxc2 = gradDensityd .* 2 .* localVxc2;
                    gradDensityVxc2{i, j, k} = localGradDensityVxc2;
                end
            end
        end

        tempVec = ElemVecToGlobal(gradDensityVxc2, HamDG.domain.numGridFine, numElem);

        Y = F * tempVec;
        ik = F.ikFine{d};
        idxnz = F.gkkFine ~= 0;
        Y(idxnz) = Y(idxnz) .* ik(idxnz);
        Y(~idxnz) = 0;
        X = real(F' * Y);

        gradGradDensityVxc2 = GlobalVecToElem(X, HamDG.domain.numGridFine, numElem);

        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    localVxc = vxc{i, j, k};
                    localGradGradDenistyVxc2 = gradGradDensityVxc2{i, j, k};
                    localVxc = localVxc - localGradGradDenistyVxc2;
                    vxc{i, j, k} = localVxc;
                end
            end
        end
    end

% -------------------- end of XC_FAMILY_GGA ----------------------------
else
    error('Unsupported XC family!');
end

Exc = Exc * HamDG.domain.Volume() / HamDG.domain.NumGridTotalFine();

end
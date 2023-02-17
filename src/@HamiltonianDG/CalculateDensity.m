function [rho, rhoLGL] = CalculateDensity(HamDG)
% HAMILTONIANDG/CALCULATEDENSITY calculate density of HamDG.
% 
%    [rho, rhoLGL] = CalculateDensity(HamDG) computes density over uniform
%    and LGL grid for each element.
%
%    See also HamiltonianDG, Spinor, SCFDG/CalculateOccupationRate.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = HamDG.numElem;
numElemTotal = prod(numElem);
grid = HamDG.grid;

occrate = HamDG.occupationRate;
numSpin = HamDG.numSpin;

rho = cell(numElemTotal, 1);
rhoLGL = cell(numElemTotal, 1);

% FIXME: there are other methods not implemented here.
% Method 3: Compute the electron density locally, and then normalize only
% in global domain, output the eigenfunctions locally.

sumRho = 0;
sumRhoLGL = 0;

% compute the local density in each element
for elemIdx = 1 : numElemTotal
    localBasis = HamDG.basisLGL{elemIdx};
    [~, numBasis] = size(localBasis);
    
    % skip the element if there is no basis functions
    if numBasis == 0
        continue;
    end
    
    % compute local wavefunction on the LGL grid
    localCoef = HamDG.eigvecCoef{elemIdx};                        
    localPsiLGL = localBasis * localCoef;
    
    % update the local density
    occrate = reshape(occrate, 1, []);
    localRhoLGL = sum(localPsiLGL.^2 .* occrate .* numSpin, 2);
                
    rhoLGL{elemIdx} = localRhoLGL;
    
    % interpolate the local density from LGL grid to the uniform grid
    transferMat = grid.LGLToUniformMatFine;
    numGridOld = grid.numLGLGridElem;
    numGridNew = grid.numUniformGridElemFine;
    localRho = InterpByTransferMat(transferMat, numGridOld, ...
                                   numGridNew, localRhoLGL);
    rho{elemIdx} = localRho;

    sumRhoLGL = sumRhoLGL + sum(localRhoLGL .* grid.LGLWeight3D(:), 'all');            
    sumRho = sumRho + sum(localRho, 'all');
end

end
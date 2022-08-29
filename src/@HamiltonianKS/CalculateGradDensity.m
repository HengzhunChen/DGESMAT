function gradDensity = CalculateGradDensity(HamKS)
% HAMILTONIANKS/CALCULATEGRADDENSITY calculates gradient of density.
% 
%    gradDensity = CalculateGradDensity(HamKS) computes gradient of 
%    density which has multi-directions for each dim and is stores as 
%    a cell structure gradDensity.
%
%    See also HamiltonianKS, HamiltonianKS/CalculateDensity.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

dim = dimDef();
gradDensity = cell(dim, 1);

F = HamKS.fft;

densityFourier = F * HamKS.density;

% compute the derivative of the density via Fourier
for d = 1 : dim
    ik = F.ikFine{d};
    idxnz = F.gkkFine ~= 0;
    
    y = densityFourier;
    y(idxnz) = y(idxnz) .* ik(idxnz);
    y(~idxnz) = 0;
    
    x = F' * y;
    gradDensity{d} = real(x);
end

end
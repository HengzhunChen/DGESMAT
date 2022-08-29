function vhart = CalculateHartree(HamKS)
% HAMILTONIANKS/CALCULATEHARTREE calculates Hartree potential.
%
%    See also HamiltonianKS, HamiltonianKS/CalculateVtot.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

F = HamKS.fft;

% The contribution of the pseudoCharge is subtracted. So the Poisson
% equation is well defined for the neutral system.

x = HamKS.density - HamKS.pseudoCharge;
y = F * x;
idxnz = F.gkkFine ~= 0;
y(idxnz) = y(idxnz) .* 4 * pi ./ F.gkkFine(idxnz);
y(~idxnz) = 0;
x = F' * y;
vhart = real(x);

end
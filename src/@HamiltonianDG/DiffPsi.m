function Dpsi = DiffPsi(HamDG, numGrid, psi, d)
% HAMILTONIANDG/DIFFPSI calculate differentiation of basis function
%
%    Dpsi = DiffPsi(HamDG, numGrid, psi, d) differentiate the basis
%    functions psi on a certain element with grids numGrid along the
%    dimension d. 
%    NOTE: numGrid should be a vector.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


DMat = HamDG.grid.DMat;

if d == 1
    vecXpsi = reshape(psi, numGrid(1), numGrid(2)*numGrid(3));
    DpsiX = DMat{1} * vecXpsi;
    Dpsi = reshape(DpsiX, numGrid(1), numGrid(2), numGrid(3));
    Dpsi = reshape(Dpsi, [], 1);
elseif d == 2
    psi = reshape(psi, numGrid(1), numGrid(2), numGrid(3));
    tempPsi = permute(psi, [2, 1, 3]);
    vecYpsi = reshape(tempPsi, numGrid(2), numGrid(1)*numGrid(3));
    DpsiY = DMat{2} * vecYpsi;
    tempPsi = reshape(DpsiY, numGrid(2), numGrid(1), numGrid(3));
    Dpsi = permute(tempPsi, [2, 1, 3]);
    Dpsi = reshape(Dpsi, [], 1);
elseif d == 3
    psi = reshape(psi, numGrid(1), numGrid(2), numGrid(3));
    tempPsi = permute(psi, [3, 2, 1]);
    vecZpsi = reshape(tempPsi, numGrid(3), numGrid(1)*numGrid(2));
    DpsiZ = DMat{3} * vecZpsi;
    tempPsi = reshape(DpsiZ, numGrid(3), numGrid(2), numGrid(1));
    Dpsi = permute(tempPsi, [3, 2, 1]);
    Dpsi = reshape(Dpsi, [], 1);
else
    error('Wrong dimension.');
end

end
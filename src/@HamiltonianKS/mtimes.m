function HX = mtimes(H, X)
% HAMILTONIANKS/MTIMES overload mutiplication operator between 
%    HamiltonianKS object H and Spinor object X.  
%
%    Usage: HX = H * X;
%    both X and HX are Spinor objects containing wave function info.
%
%    See also HamiltonianKS, Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


F = H.fft;

if ~F.isInitialized || ~F.isInitializedFine
    error('Fourier is not prepared.');
end


ntot = H.domain.NumGridTotal();
ntotFine = H.domain.NumGridTotalFine();
vol = H.domain.Volume();

XFine = CoarseToFine(X, F);

% add local pseudopotential
HX_p = XFine .* H.vtot;

% add nonlocal pseduopotential
pseudoList = H.pseudoList;
natom = length(pseudoList);
for i = 1 : natom
    vnlList = pseudoList(i).vnlList;
    nobt = length(vnlList);
    for j = 1 : nobt
        vnl = vnlList(j);
        vnlwgt = vnl.wgt;
        iv = vnl.idx;
        dv = vnl.val;

        weight = dv(:, 1)' * XFine(iv, :);  % index 1 means value component
        weight = weight * vol / ntotFine * vnlwgt;
        
        HX_p(iv, :) = HX_p(iv, :) + dv(:, 1) .* weight;
    end
end

% add kinetic part in Fourier space and then back to real space
HX_kg = F * X .* F.gkk ./ 2;

fac = sqrt( ntotFine / ntot );
idx = F.idxFineGrid;
yfine = F * HX_p;
x = yfine(idx, :) .* fac;

HXg = HX_kg + x;
HX = F' * HXg;
HX = real(HX);  % important


if H.isHybrid && H.isEXXActive
    % TODO
end


% apply filter on the wavefunctions before exit, if required
if H.ecutFilter.isApplyFilter == 1
    HY = F * HX;
    idx = (F.gkk ./ 2) > H.ecutFilter.wfnCutoff;
    HY(idx, :) = 0;
    HX = F' * HY;
    HX = real(HX);
end

end
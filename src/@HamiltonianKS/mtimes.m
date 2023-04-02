function HX = mtimes(H, X)
% HAMILTONIANKS/MTIMES overload mutiplication operator between 
%    HamiltonianKS object H and Spinor object X.  
%
%    Usage: HX = H * X;
%    both X and HX are Spinor objects containing wave function info.
%
%    See also HamiltonianKS, Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


F = H.fft;
if ~F.isInitialized || ~F.isInitializedFine
    error('Fourier is not prepared.');
end

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
    if ~isempty(vnlList)
        iv = vnlList.idx;
        vnl = vnlList.val;
        wgt = vnlList.wgt .* (vol / ntotFine);
        wgt = wgt .* (vnl' * XFine(iv, :));
        HX_p(iv, :) = HX_p(iv, :) + vnl * wgt;
    end
end

% add kinetic part in Fourier space and then back to real space
HX_kg = F * X .* F.gkk ./ 2;

idx = F.idxFineGrid;
yfine = F * HX_p;
x = yfine(idx, :);

HXg = HX_kg + x;
HX = F' * HXg;
HX = real(HX);  % important


% apply filter on the wavefunctions before exit, if required
if H.ecutFilter.isApplyFilter == 1
    HY = F * HX;
    idx = (F.gkk ./ 2) > H.ecutFilter.wfnCutoff;
    HY(idx, :) = 0;
    HX = F' * HY;
    HX = real(HX);
end

end
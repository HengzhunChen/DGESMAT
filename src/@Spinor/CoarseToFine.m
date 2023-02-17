function Xfine = CoarseToFine(X, F)
% SPINOR/COARSETOFINE map wave function from coarse grid to fine grid 
%    using Fourier object F.
%
%    Xfine = CoarseToFine(X, F) returns the Spinor object Xfine over 
%    fine grid by using the Fourier transform of wave function over 
%    coarse grid.
%
%    See also Spinor, Fourier, Domain.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

ntot = F.domain.NumGridTotal();
ntotFine = F.domain.NumGridTotalFine();

x = X.wavefun;
[~, ncols] = size(x);
yfine = zeros(ntotFine, ncols);

fac = sqrt( ntot / ntotFine );
idx = F.idxFineGrid;

y = F * x;
yfine(idx, :) = y .* fac;
xfine = F' * yfine;

Xfine = Spinor(X.domain, X.numStateList, real(xfine));

end
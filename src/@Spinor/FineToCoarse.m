function X = FineToCoarse(Xfine, F)
% SPINOR/FINETOCOARSE map wave function from fine grid to coarse grid 
%    using Fourier object F.
%
%    X = FineToCoarse(Xfine, F) returns the Spinor object X over 
%    coarse grid by using the Fourier transform of wave function over 
%    fine grid.
%
%    See also Spinor, Fourier, Domain.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

idx = F.idxFineGrid;
yfine = F * Xfine;
X = F' * yfine(idx, :);
X = real(X);

end
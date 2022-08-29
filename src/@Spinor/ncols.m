function w = ncols(X)
% SPINOR/NCOLS Number of columns of the wave functions
%
%    w = ncols(X) returns the number of columns in wavefun.
%
%    See also Spinor, nrows.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

w = size(X.wavefun, 2);

end
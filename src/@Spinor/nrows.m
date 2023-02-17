function h = nrows(X)
% SPINOR/NROWS Number of rows of the wave functions
%
%    h = nrows(X) returns the number of rows in wavefun.
%
%    See also Spinor, ncols.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

h = size(X.wavefun, 1);

end
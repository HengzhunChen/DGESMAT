function X = real(X)
% SPINOR/REAL Real function for wave functions
%
%    X = real(X) returns a Spinor object with wave function of the real 
%    part of X.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

X.wavefun = real(X.wavefun);

end
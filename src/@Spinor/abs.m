function Xabs = abs(X)
% SPINOR/ABS Abs function for wave function.
%
%    Xabs = abs(X) returns a Spinor object with the absolute value of the 
%    wave function of X.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

Xabs = X;
Xabs.wavefun = abs(X.wavefun);

end
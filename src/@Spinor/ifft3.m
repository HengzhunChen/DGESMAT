function X = ifft3(Y)
% SPINOR/FFT3 ifft3 function for wave function
%
%    X = ifft3(Y) returns a Spinor object with the inverse fast Fourier 
%    transform of the wave function of X.
%
%    See also Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

numGrid = Y.domain.numGrid;
vol = Y.domain.Volume();

y = Y.wavefun;
x = ifft3Execute(y, numGrid, vol);
X = Y;
X.wavefun = x;

end
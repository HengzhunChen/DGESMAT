function Y = fft3(X)
% SPINOR/FFT3 fft3 function for wave function
%
%    Y = fft3(X) returns a Spinor object with the fast Fourier transform 
%    of the wave function of X.
%
%    See also Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

numGrid = X.domain.numGrid;
vol = X.domain.Volume();

x = X.wavefun;
y = fft3Execute(x, numGrid, vol);
Y = X;
Y.wavefun = y;

end
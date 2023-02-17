function y = fft3Execute(x, numGrid, vol)
% FFT3EXECUTE Implementation of 3-dim fast Fourier transformation
%
%    y = fft3Execute(x, numGrid, vol) returns the fast Fourier tranform of
%    x whose is 3-dim tensor (stores as matrix numGridTotal * numCol) over 
%    domain with number of grid numGrid (vector) and volume vol.
%
%    See also ifft3Execute, Fourier, Spinor/fft3.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

ntot = prod(numGrid);
factor = vol / ntot;

% method 1
ncols = size(x, 2);
X = reshape(x, [numGrid, ncols]);
X = fft(X, [], 1);
X = fft(X, [], 2);
X = fft(X, [], 3);
y = reshape(X, [ntot, ncols]);

% method 2
% [nrows, ncols] = size(x);
% y = zeros(nrows, ncols);
% for j = 1 : ncols
%     X = reshape(x(:, j), numGrid);
%     Y = fftn(X);
%     y(:, j) = reshape(Y, [ntot, 1]);
% end

y = y .* factor;

end
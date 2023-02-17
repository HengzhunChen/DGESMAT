function y = ifft3Execute(x, numGrid, vol)
% IFFT3EXECUTE Implementation of 3-dim inverse fast Fourier transformation
%
%    y = ifft3Execute(x, numGrid, vol) returns the inverse fast Fourier 
%    tranform of x whose is 3-dim tensor (stores as matrix 
%    numGridTotal * numCol) over domain with number of grid numGrid 
%    (vector) and volume vol.
%
%    See also fft3Execute, Fourier, Spinor/ifft3.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

ntot = prod(numGrid);
ifactor = ntot / vol;

% method 1
ncols = size(x, 2);
X = reshape(x, [numGrid, ncols]);
X = ifft(X, [], 1);
X = ifft(X, [], 2);
X = ifft(X, [], 3);
y = reshape(X, [ntot, ncols]);

% method 2
% [nrows, ncols] = size(x);
% y = zeros(nrows, ncols);
% for j = 1 : ncols
%     X = reshape(x(:, j), numGrid);
%     Y = ifftn(X);
%     y(:, j) = reshape(Y, [ntot, 1]);
% end

y = y .* ifactor;

end
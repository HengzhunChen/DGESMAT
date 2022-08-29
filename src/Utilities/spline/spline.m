function [b, c, d] = spline(x, y)
% SPLINE cubic interpolating spline
%
%    [b, c, d] = spline(x, y) computes the cubic interpolation spline
%    coefficients b, c, d according to the arrays of data abscissas and 
%    ordinates x, y.
%
%    This SPLINE function is designed specifically for the interpolation
%    part for pseudopotential generation in the electronic structure
%    calculation.
%    The coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
%    for a cubic interpolating spline
%
%     s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
%
%    for  x(i) .le. x .le. x(i+1)
%
%    NOTE:
%    (1) x is the abscissas of the knots and must in a strictly increasing 
%        order.
%    (2) This cubic interpolation spline uses end condition utilizing the 
%        third derivatives.
% 
%    See also seval, spline_rad.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


n = length(x);
x = reshape(x, n, 1);
y = reshape(y, [], 1);

b = zeros(n, 1);
c = zeros(n, 1);
d = zeros(n, 1);

if n < 2
    error(" SPLINE requires n >= 2!");
elseif n == 2
    b(1) = (y(2) - y(1)) / (x(2) - x(1));
    b(2) = b(1);
    return;
end

% set up tridiagonal system
% b = diagonal, d = offdiagonal, c = right hand side
d(1 : n-1) = x(2 : n) - x(1 : n-1);
b(2 : n-1) = 2 * (d(1 : n-2) + d(2 : n-1));
dy = y(2 : n) - y(1 : n-1);
c(2 : n-1) = dy(2 : n-1) ./ d(2 : n-1) - dy(1 : n-2) ./ d(1 : n-2);

% end conditions. Third derivatives at x(1) and x(n) obtained from divided
% differences.
b(1) = -d(1);
b(n) = -d(n-1);
c(1) = 0;
c(n) = 0;
if n > 3
    c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1));
    c(1) = c(1) * d(1) * d(1) / (x(4)-x(1));
    c(n) = c(n-1)/(x(n)- x(n-2)) - c(n-2)/(x(n-1)-x(n-3));
    c(n) = -c(n) * d(n-1) * d(n-1) / (x(n) - x(n-3));
end

% forward elimination
for i = 2 : n
    t = d(i-1) / b(i-1);
    b(i) = b(i) - t * d(i-1);
    c(i) = c(i) - t * c(i-1);
end

% backward substitution
c(n) = c(n) / b(n);
for i = n-1 : -1 : 1
    c(i) = (c(i) - d(i) * c(i+1)) / b(i);
end

% compute polynomial coefficients
b(n) = (y(n) - y(n-1)) / d(n-1) + d(n-1) * (c(n-1) + 2*c(n));

b(1 : n-1) = (y(2 : n) - y(1 : n-1)) ./ d(1 : n-1) ...
                - d(1 : n-1) .* ( c(2 : n) + 2 .* c(1 : n-1) );
d(1 : n-1) = (c(2 : n) - c(1 : n-1)) ./ d(1 : n-1);
c(1 : n-1) = 3 * c(1 : n-1);

c(n) = 3 * c(n);
d(n) = d(n-1);

end
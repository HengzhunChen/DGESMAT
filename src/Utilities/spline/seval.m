function v = seval(u, x, y, b, c, d)
% SEVAL evaluates the spline function
%
%    v = seval(u, x, y, b, c, d) evaluates the spline function value v at 
%    the abscissa u according to spline data x, y, b, c, d.
%    x, y = the arrays of data abscissas and ordinates,
%    b, c, d = arrays of spline coefficients computed by spline.
%
%    This SEVAL function is designed specifically for the interpolation
%    part for pseudopotential generation in the electronic structure
%    calculation. Therefore if u is outside the range [min(x), max(x)],
%    the corresponding v value will be an extrapolation.
%    It evaluates the spline function acccording to 
% 
%     seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
% 
%    where x(i) .lt. u .lt. x(i+1), using horner's rule
% 
%     if  u .lt. x(1) then  i = 1  is used.
%     if  u .ge. x(n) then  i = n  is used.
%
%     See also spline.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


n = length(x);
if n < 2
    error(' spline requires n >= 2!');
end

x = reshape(x, [], 1);  % make sure x, y are column vectors
y = reshape(y, [], 1);  

breaks = [x; Inf];
coefs = [d, c, b, y];
pp = mkpp(breaks, coefs);
v = ppval(u, pp);


end
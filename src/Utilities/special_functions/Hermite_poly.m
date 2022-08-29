function y = Hermite_poly(x, order)
% HERMITE_POLY evaluate Hermite polynomial
%
%    y = Hermite_poly(x, order) evaluates order-th order Hermite 
%    polynomials (physicist convention) value y at x.
%

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

switch order
    case 0
        y = 1.0 * ones(size(x));
    case 1
        y = 2.0 * x;
    case 2
        y = 4.0 * x.^2 - 2.0;
    case 3 
        y = 8.0 * x.^3 - 12.0 * x;
    case 4 
        y = 16.0 * x.^4 - 48.0 * x.^2 + 12.0;
    case 5 
        y = 32.0 * x.^5 - 160.0 * x.^3 + 120.0 * x;
    case 6 
        y = 64.0 * x.^6 - 480.0 * x.^4 + 720.0 * x.^2 - 120.0;
end

end
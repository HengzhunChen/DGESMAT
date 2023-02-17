function [epsc, vc] = Correlation_pz(rs)
% CORRELATION_PZ subroutine to calculate PZ type correlation
%
%    See also xcRef.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

a = 0.0311;
b = -0.048;
c = 0.0020;
d = -0.0116;
gc = -0.1423;
b1 = 1.0529;
b2 = 0.3334;

vc = zeros(size(rs));
epsc = zeros(size(rs));

% high density formula
idxl = rs < 1;
rsl = rs(idxl);
lnrs = log (rsl);
epsc(idxl) = a*lnrs + b + c*rsl.*lnrs + d*rsl;
vc(idxl) = a*lnrs + (b-a/3) + 2/3*c*rsl.*lnrs + (2*d-c)/3*rsl;

% interpolation formula
idxg = rs >= 1;
rsg = rs(idxg);
rs12 = sqrt(rsg);
ox = 1 + b1*rs12 + b2*rsg;
dox = 1 + 7/6*b1*rs12 + 4/3*b2*rsg;
epsc(idxg) = gc./ ox;
vc(idxg) = epsc(idxg) .* dox ./ ox;

end
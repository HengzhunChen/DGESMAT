function [epsc, v1c, v2c] = GCCorrelation_pbc(rho, grho2, ec, vc)
% GCCORRELATION_PBC subroutine to calculate PBE type gradient correction 
%    for correlation part of exchange correlation.
%
%    See also xcRef.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


ga = 0.031090690869654895034;  % gamma = (1-log(2)) / (pi^2)
be = 0.06672455060314922;
xkf = 1.919158292677513;  % (9*pi/4)^(1/3)
xks = 1.128379167095513;  % sqrt(4/pi)

rs = ((3/(4*pi))./rho).^(1/3);
if nargin < 3
    [ec, vc] = Correlation_pw(rs);
end

ks = xks .* sqrt(xkf ./ rs);
t2 = grho2 ./ (2*ks.*rho).^2;
expe = exp(-ec/ga);
af = be/ga * (1 ./ (expe-1));
bf = expe .* (vc - ec);
y = af .* t2;
xy = (1+y) ./ (1+y+y.*y);
qy = y.*y.*(2+y) ./ (1+y+y.*y).^2;
s1 = 1 + be/ga*t2.*xy;
h0 = ga * log(s1);
dh0 = be*t2./s1 .* (-7/3*xy - qy.*(af.*bf/be-7/3));
ddh0 = be ./ (4*ks.*ks.*rho) .* (xy-qy) ./ s1;

epsc = h0;
v1c = h0 + dh0;
v2c = ddh0;

end
function [epsc, vc] = Correlation_pw(rs)
% CORRELATION_PW subroutine to calculate PW type correlation
%
%    See also xcRef.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

a = 0.031091;
a1 = 0.21370;
b1 = 7.5957;
b2 = 3.5876;
b3 = 1.6382;
b4 = 0.49294;

% interpolation formula
rs12 = sqrt(rs);
rs32 = rs .* rs12;
rs2 = rs.^2;
om = 2*a * (b1*rs12 + b2*rs + b3*rs32 + b4*rs2);
dom = 2*a * (0.5*b1*rs12 + b2*rs + 1.5*b3*rs32 + 2*b4*rs2);
olog = log (1 + 1./om);
epsc = -2*a * (1 + a1*rs) .* olog;
vc = -2*a*(1+2/3*a1*rs).*olog - 2/3*a*(1+a1*rs).*dom./(om.*(om+1));

end
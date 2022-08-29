function [epsc, v1c, v2c] = GGA_C_PBE(rho, grho2)
% GGA_C_PBE calculates correlation part of GGA type in PBE functionals
%
%    See also xcRef, xc_gga_exc_vxc, GGA_X_PBE.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

% NOTE:
% pbe = sla + pw + pbx + pbc
% pbe_c = pw + pbc

v1c = zeros(size(rho));
v2c = zeros(size(rho));
epsc = zeros(size(rho));

idxxc = abs(rho) > 1e-10;
rhoxc = abs(rho(idxxc));

rs = ((3/(4*pi))./rhoxc).^(1/3);
[epsc_pw, vc_pw] = Correlation_pw(rs);

v1c(idxxc) = vc_pw;
epsc(idxxc) = epsc_pw;

idxcxc = abs(rho) > 1e-6 & grho2 > 1e-10;
rhocxc = abs(rho(idxcxc));
grho2cxc = grho2(idxcxc);

[epsc_gc, v1c_gc, v2c_gc] = GCCorrelation_pbc(rhocxc, grho2cxc);

v1c(idxcxc) = v1c(idxcxc) + v1c_gc;
epsc(idxcxc) = epsc(idxcxc) + epsc_gc;

v2c(idxcxc) = v2c_gc;

end
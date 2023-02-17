function [epsx, v1x, v2x] = GGA_X_PBE(rho, grho2)
% GGA_X_PBE calculates exchange part of GGA type in PBE functionals
%
%    See also xcRef, xc_gga_exc_vxc, GGA_C_PBE.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

% NOTE:
% pbe = sla + pw + pbx + pbc
% pbe_x = sla + pbx

v1x = zeros(size(rho));
v2x = zeros(size(rho));
epsx = zeros(size(rho));

idxxc = abs(rho) > 1e-10;
rhoxc = abs(rho(idxxc));

rs = ((3/(4*pi))./rhoxc).^(1/3);
[epsx_sla, vx_sla] = Exchange_sla(rs);

v1x(idxxc) = vx_sla;
epsx(idxxc) = epsx_sla;

idxcxc = abs(rho) > 1e-6 & grho2 > 1e-10;
rhocxc = abs(rho(idxcxc));
grho2cxc = grho2(idxcxc);

[epsx_gc, v1x_gc, v2x_gc] = GCExchange_pbx(rhocxc, grho2cxc);

v1x(idxcxc) = v1x(idxcxc) + v1x_gc;
epsx(idxcxc) = epsx(idxcxc) + epsx_gc;

v2x(idxcxc) = v2x_gc;

end
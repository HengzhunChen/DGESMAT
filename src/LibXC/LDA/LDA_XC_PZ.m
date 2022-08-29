function [epsxc, vxc] = LDA_XC_PZ(rho)
% LDA_XC_PZ computes LDA type exchange correlation in PZ functional.
%
%    See also xcRef, xc_lda_exc_vxc.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

vxc = zeros(size(rho));
epsxc = zeros(size(rho));

idx = rho > 1e-10;
rhoxc = rho(idx);

rs = ((3/(4*pi)) ./ rhoxc) .^ (1/3);
[epsx, vx] = Exchange_sla(rs);
[epsc, vc] = Correlation_pz(rs);

vxc(idx) = vx + vc;
epsxc(idx) = epsx + epsc;

end
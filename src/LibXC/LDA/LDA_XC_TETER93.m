function [epsxc, vxc] = LDA_XC_TETER93(rho)
% LDA_XC_TETER93 computes LDA type exchange correlation in TETER93
% functional.
%
%    See also xcRef, xc_lda_exc_vxc.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


a0 = 0.4581652932831429;
a1 = 2.217058676663745;
a2 = 0.7405551735357053;
a3 = 0.01968227878617998;
b1 = 1.0;
b2 = 4.504130959426697;
b3 = 1.110667363742916;
b4 = 0.02359291751427506;

rs = ( 3/(4*pi) ./ rho ) .^ (1/3);

frac_up = a0 + rs .* ( a1 + rs .* (a2 + rs .* a3) );
frac_down = rs .* ( b1 + rs .* (b2 + rs .* (b3 + rs .* b4)) );

epsxc = - frac_up ./ frac_down;

dExcdr = -( (a1 + 2*a2.*rs + 3*a3.*rs.^2) .* frac_down - ... 
        frac_up .* (b1 + 2*b2.*rs + 3*b3.*rs.^2 + 4*b4.*rs.^3) ) ...
        ./ (frac_down.^2);
vxc = epsxc - (1/3) .* dExcdr .* rs; 


end
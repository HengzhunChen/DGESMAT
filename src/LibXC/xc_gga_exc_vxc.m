function [epsxc, vrho, vsigma] = xc_gga_exc_vxc(XCFuncType, rho, grho2)
% XC_LDA_EXC_VXC Interface of different types LDA exchange correlation
%
%    See also xcRef.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% ======================================================================
% GGA exchange-correlation
% ======================================================================

if XCFuncType == "XC_GGA_XC_PBE"
    [epsx, vrhox, vsigmax] = GGA_X_PBE(rho, grho2);
    [epsc, vrhoc, vsigmac] = GGA_C_PBE(rho, grho2);    
    epsxc = epsx + epsc;
    vrho = vrhox + vrhoc;
    vsigma = vsigmax + vsigmac;
    return;
else
    % TODO
end


% ======================================================================
% hybrid GGA functionals
% ======================================================================

if XCFuncType == "XC_HYB_GGA_XC_HSE06"
    [epsxc, vrho, vsigma] = HYB_GGA_XC_HSE06(rho, grho2);
    return;
else
    error('Unsupported XC type');
end


end
function [epsxc, vxc] = xc_lda_exc_vxc(XCFuncType, rho)
% XC_LDA_EXC_VXC Interface of different types LDA exchange correlation
%
%    See also xcRef.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% ======================================================================
% LDA exchange-correlation
% ======================================================================

if XCFuncType == "XC_LDA_XC_TETER93"
    [epsxc, vxc] = LDA_XC_TETER93(rho);
elseif XCFuncType == "XC_LDA_XC_PZ"
    [epsxc, vxc] = LDA_XC_PZ(rho);
else
    error('Unsupported XC Type');
end


% ======================================================================
% hybrid LDA functionals
% ======================================================================

% TODO


end
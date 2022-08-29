function vtot = CalculateVtot(HamKS)
% HAMILTONIANKS/CALCULATEVTOT calculates total local potential
%
%    See also HamiltonianKS.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

global esdfParam;

if ~esdfParam.userOption.general.isUseVLocal
    vtot = HamKS.vext + HamKS.vhart + HamKS.vxc;
else
    vtot = HamKS.vext + HamKS.vLocalSR + HamKS.vhart + HamKS.vxc;
end

end
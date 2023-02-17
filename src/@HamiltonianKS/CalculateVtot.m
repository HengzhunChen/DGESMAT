function vtot = CalculateVtot(HamKS)
% HAMILTONIANKS/CALCULATEVTOT calculates total local potential
%
%    See also HamiltonianKS.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


if ~HamKS.userOption.isUseVLocal
    vtot = HamKS.vext + HamKS.vhart + HamKS.vxc;
else
    vtot = HamKS.vext + HamKS.vLocalSR + HamKS.vhart + HamKS.vxc;
end

end
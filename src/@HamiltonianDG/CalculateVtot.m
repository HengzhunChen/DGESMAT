function vtot = CalculateVtot(HamDG)
% HAMILTONIANDG/CALCULATEVTOT calculates total local potential
%
%    vtot = CalculateVtot(HamDG) returns a cell structure containing total
%    local potential for each element.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

numElemTotal = prod(HamDG.numElem);
vtot = cell(numElemTotal, 1);

for elemIdx = 1 : numElemTotal
    localVext = HamDG.vext{elemIdx};
    localVhart = HamDG.vhart{elemIdx};
    localVxc = HamDG.vxc{elemIdx};
    
    vtot{elemIdx} = localVext + localVhart + localVxc;
end

end
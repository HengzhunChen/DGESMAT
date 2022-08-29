function vtot = CalculateVtot(HamDG)
% HAMILTONIANDG/CALCULATEVTOT calculates total local potential
%
%    vtot = CalculateVtot(HamDG) returns a cell structure containing total
%    local potential for each element.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

vtot = cell(HamDG.numElem);

for k = 1 : HamDG.numElem(3)
    for j = 1 : HamDG.numElem(2)
        for i = 1 : HamDG.numElem(1)
            localVext = HamDG.vext{i, j, k};
            localVhart = HamDG.vhart{i, j, k};
            localVxc = HamDG.vxc{i, j, k};
            
            vtot{i, j, k} = localVext + localVhart + localVxc;
        end
    end
end

end
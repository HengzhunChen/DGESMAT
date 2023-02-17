function maxForce = MaxForce(atomList)
% MAXFORCE computes the maximum of magnitude of force in atomList 
%
%    See also Atom.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

numAtom = length(atomList);
maxForce = 0;

for i = 1 : numAtom
    forceMag = norm(atomList(i).force, 2);
    maxForce = max(maxForce, forceMag);
end

end
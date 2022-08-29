function gridpos = UniformMesh(domain)
% UNIFORMMESH Generate a uniform mesh from a domain
%
%    gridpos = UniformMesh(domain) returns a uniform grid mesh gridpos
%    over domain.
%    NOTE: gridpos is a row vectors
%
%    See also Domain, UniformMeshFine.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

dim = dimDef();
gridpos = cell(dim, 1);
for d = 1 : dim
    h = domain.length(d) / domain.numGrid(d);
    gridpos{d} = domain.posStart(d) + h * (0 : domain.numGrid(d)-1);
end

end
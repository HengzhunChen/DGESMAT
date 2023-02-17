function  gridpos = UniformMeshFine(domain)
% UNIFORMMESHFINE Generate a uniform fine mesh from a domain
%
%    gridpos = UniformMeshFine(domain) returns a fine uniform grid mesh 
%    gridpos over domain.
%    NOTE: gridpos is a row vectors.
%
%    See also Domain, UniformMesh.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

dim = dimDef();
gridpos = cell(dim, 1);
for d = 1 : dim
    h = domain.length(d) / domain.numGridFine(d);
    gridpos{d} = domain.posStart(d) + h * (0 : domain.numGridFine(d)-1);
end

end

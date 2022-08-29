function gridpos = LGLMesh(domain, numLGLGrid)
% LGLMESH Generate a LGL mesh from a domain.
%
%    gridpos = LGLMesh(domain, numLGLGrid) returns a LGL grid mesh over
%    domain according to number of LGL grid numLGLGrid(vector).
%    NOTE: gridpos is a row vectors.
%
%    See also Domain, GenerateLGL.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


dim = dimDef();
gridpos = cell(dim, 1);

for d = 1 : dim
    [mesh, ~, ~, ~] = GenerateLGL( numLGLGrid(d) );
    mesh = mesh';
    gridpos{d} = domain.posStart(d) + (mesh + 1.0) * domain.length(d) * 0.5;
end

end
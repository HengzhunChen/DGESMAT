function atomDensity = CalculateAtomDensity(PT, atom, domain, gridpos)
% PERIODTABLE/CALCULATEATOMDENSITY generates the atomic density. 
%    This is to be used with structure factor to generate the initial 
%    atomic density, which can be used both for PWDFT and DGDFT.
%
%    See also PeriodTable, HamiltonianKS/CalculateAtomDensity,
%    HamiltonianDG/CalculateAtomDensity.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% NOTE: This assumes atom.pos should agree with posStart in domain
% (PWDFT) and global domain (DGDFT).

dim = dimDef();
L = domain.length;

spldata = PT.splineMap(atom.type);
Rzero = PT.RcutRhoAtom(atom.type);

numGrids = length(gridpos{1}) * length(gridpos{2}) * length(gridpos{3});
atomDensity = zeros(numGrids, 1);

% Compute the minimal distance of the atom to this set of grid points and
% determine whether to continue (at least one grid within Rzero)
dist = cell(dim, 1);
minDist = zeros(dim, 1);
for d = 1 : dim
    dist{d} = gridpos{d} - atom.pos(d);
    dist{d} = dist{d} - round( dist{d} / L(d) ) * L(d);
    minDist(d) = min( [abs(dist{d}), Rzero] );
end

if norm(minDist) <= Rzero
    % At least one grid point is within Rzero
    
    % compute distance between grid points and atom
    [distI, distJ, distK] = ndgrid(dist{1}, dist{2}, dist{3});
    distxyz = [distI(:), distJ(:), distK(:)];
    distTemp = sqrt( sum(distxyz .* distxyz, 2) );
    
    idx = find( distTemp <= Rzero );
    rad = distTemp(idx);
    
    idxsize = length(idx);
    if idxsize > 0
        valspl = spldata.rhoAtom;
        val = seval(rad, valspl(:, 1), valspl(:, 2), valspl(:, 3), valspl(:, 4), valspl(:, 5));
        atomDensity(idx) = val; 
    end
end

end
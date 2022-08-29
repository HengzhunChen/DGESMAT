function pseudoCharge = CalculatePseudoCharge(PT, atom, domain, gridpos)
% PERIODTABLE/CALCULATEPSEUDOCHARGE generates the pseudo-charge and its 
%    derivatives.
%     
%    pseudoCharge = CalculatePseudoCharge(PT, atom, domain, gridpos)
%    returns a struct pseudoCharge whose members are idx, val (see also
%    NOTE below). To Check whether pseudoCharge is empty, use
%    `isempty(pseudoCharge.idx)` instead of `isempty(pseudoCharge)`.
%    This function is used WHEN user option isUseVLocal is false.
%
%    See also PeriodTable, HamiltonianKS/CalculatePseudoPotential,
%    HamiltonianDG/CalculatePseudoPotential.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% The start of the radial grid
MIN_RADIAL = 1e-10;

dim = dimDef();
L = domain.length;

spldata = PT.splineMap(atom.type);
Rzero = PT.RcutPseudoCharge(atom.type);

pseudoCharge = struct(...
    'idx', [], ...  % index within cutoff
    'val', [] ...   % value and its three derivatives
    );
% NOTE: to check whether pseudoCharge is empty, use 
% isempty(pseudoCharge.idx) instead of isempty(pseudoCharge)


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
    xx = distxyz(idx, 1);
    yy = distxyz(idx, 2);
    zz = distxyz(idx, 3);
    
    idxsize = length(idx);
    
    valspl = spldata.pseudoCharge;
    val = seval(rad, valspl(:, 1), valspl(:, 2), valspl(:, 3), valspl(:, 4), valspl(:, 5));
    
    derspl = spldata.drv_pseudoCharge;
    der = seval(rad, derspl(:, 1), derspl(:, 2), derspl(:, 3), derspl(:, 4), derspl(:, 5));
    
    dv = zeros(idxsize, dim+1);  % value and its three derivatives
    idxnz = rad > MIN_RADIAL;
    dv(:, 1) = val;  % value
    dv(idxnz, 2) = der(idxnz) .* xx(idxnz) ./ rad(idxnz);  % DX
    dv(idxnz, 3) = der(idxnz) .* yy(idxnz) ./ rad(idxnz);  % DY
    dv(idxnz, 4) = der(idxnz) .* zz(idxnz) ./ rad(idxnz);  % DZ

    pseudoCharge.idx = idx;
    pseudoCharge.val = dv;
end

end
function [VLocalSR, GaussianPseudoCharge] = CalculateVLocal(PT, atom, domain, gridpos)
% PERIODTABLE/CALCULATEVLOCAL generates the short range local pseudo 
%    potential and its derivatives, as well as the Gaussian compensation 
%    charge and its derivatives.
% 
%    [VLocalSR, GaussianPseudoCharge] = CalculateVLocal(PT, atom, domain, gridpos)
%    Both VLocalSR and GaussianPseudoCharge are struct with idx and val. 
%    To check whether they are empty or not, use isempty to check it
%    member like idx instead of itself.
%    This function is used WHEN user option isUseVLocal is true.
%
%    See also PeriodTable, PeriodTable/Setup, 
%    HamiltonianKS/CalculatePseudoPotential, 
%    HamiltonianDG/CalculatePseudoPotential.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


dim = dimDef();
L = domain.length;  % vector

spldata = PT.splineMap(atom.type);

% use the pseudocharge cutoff for Gaussian compensation charge and short
% range potential
Rzero = PT.RcutPseudoCharge(atom.type);
RGaussian = PT.RGaussian(atom.type);
Zion = PT.Zion(atom.type);


VLocalSR = struct(...
    'idx', [], ...  % index within cutoff
    'val', [] ...   % value and its three derivatives
    );

GaussianPseudoCharge = struct(...
    'idx', [], ...  % index within cutoff
    'val', [] ...   % value and its three derivatives
    );
% NOTE: to check whether the above struct is empty, use 
% isempty(VLocalSR.idx) instead of isempty(VLocal)


% compute the minimal distance of the atom to this set of grid points and
% determine whether to continue
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
    
    % -------------- short range pseudopotential ----------------------
    valspl = spldata.vLocal;
    val = seval(rad, valspl(:, 1), valspl(:, 2), valspl(:, 3), valspl(:, 4), valspl(:, 5));
    
    derspl = spldata.drv_vLocal;
    der = seval(rad, derspl(:, 1), derspl(:, 2), derspl(:, 3), derspl(:, 4), derspl(:, 5));
    
    dv = zeros(idxsize, dim+1);  % value and its three derivatives
    dv(:, 1) = val;
    % FIXME: derivatives later
    
    VLocalSR.idx = idx;
    VLocalSR.val = dv;
    
    % ------------------ Gaussian pseudocharge ------------------------
    factor = Zion / ( sqrt(pi) * RGaussian )^3;
    dv = zeros(idxsize, dim+1);
    dv(:, 1) = factor * exp( -(rad / RGaussian).^2 );
    % FIXME: derivatives later
    
    GaussianPseudoCharge.idx = idx;
    GaussianPseudoCharge.val = dv;
    
end
    
end
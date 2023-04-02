function vnlList = CalculateNonlocalPP(PT, atom, domain, gridpos)
% PERIODTABLE/CALCULATENONLOCALPP generates the nonlocal pseudopotential 
%     projectors.
%
%    vnlList = CalculateNonlocalPP(PT, atom, domain, gridpos) returns a
%    struct vnlList whose members include idx, val, drv and wgt. idx is
%    index of nonlocal projectors, val and drv are value and three
%    derivatives of nonlocal projectors, wgt is weight of nonlocal
%    projectors.
%
%    See also PeriodTable, HamiltonianKS/CalculatePseudoPotential,
%    HamiltonianDG/CalculatePseudoPotential.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% TODO: 
% move spherical harmonic function into utility to simplify the code

% The start of the radial grid
MIN_RADIAL = 1e-10;


dim = dimDef();
L = domain.length;  % vector

% get entry data and spline data
ptentry = PT.pteMap(atom.type);
spldata = PT.splineMap(atom.type).nonlocal;

Rzero = PT.RcutNonlocal(atom.type);

% Initialize 
% first count all the pseudopotentials
numpp = 0;
n = length(ptentry.nltypes);
for g = 1 : n
    typ = ptentry.nltypes(g);
    if typ == 0
        numpp = numpp + 1;
    elseif typ == 1
        numpp = numpp + 3;
    elseif typ == 2
        numpp = numpp + 5;
    elseif typ == 3
        numpp = numpp + 7;
    end
end

vnlList = struct(...
    'idx', [], ...  % index within cutoff
    'val', [], ...  % value of nonlocal potential projector
    'drv', [], ...  % derivatives in 3 directions
    'wgt', [] ...  % weight of nonlocal projectors
);
% NOTE: in current implementation, all nonlocal orbits share the same
% nonzero index within cutoff, see also PeriodTable/RcutNonlocal.


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
    vnlList.idx = idx;
    vnlList.val = zeros(idxsize, numpp);
    vnlList.wgt = zeros(numpp, 1);
    vnlList.drv = cell(numpp, 1);

    % process non-local pseudopotential one by one
    countpp = 0;
    n = length(ptentry.nltypes);
    for g = 1 : n
        wgt = ptentry.nlweights(g);
        typ = ptentry.nltypes(g);
        
        idxStart = 10*(g-1) + 1;

        valspl = spldata(:, idxStart:idxStart+4);
        val = seval(rad, valspl(:, 1), valspl(:, 2), valspl(:, 3), valspl(:, 4), valspl(:, 5));
    
        derspl = spldata(:, idxStart+5 : idxStart+9);
        der = seval(rad, derspl(:, 1), derspl(:, 2), derspl(:, 3), derspl(:, 4), derspl(:, 5));

        if typ == PT.ptNLtype.L0
            coef = sqrt(1 / (4 * pi));  % spherical harmonics
            dv = zeros(idxsize, dim);  % value of three derivatives
            
            idxnz = rad > MIN_RADIAL;
            vnlval = coef * val;  % value
            dv(idxnz, 1) = coef * der(idxnz) .* xx(idxnz) ./ rad(idxnz);  % DX
            dv(idxnz, 2) = coef * der(idxnz) .* yy(idxnz) ./ rad(idxnz);  % DY
            dv(idxnz, 3) = coef * der(idxnz) .* zz(idxnz) ./ rad(idxnz);  % DZ
            
            countpp = countpp + 1;
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
        end
        
        if typ == PT.ptNLtype.L1
            coef = sqrt(3 / (4 * pi));  % spherical hamonics

            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors            
            dv = zeros(idxsize, dim);  % value of three derivatives
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            vnlval(k) = coef *( (xx(k) ./ rad(k)) .* val(k) );
            dv(k, 1) = coef *( (der(k) - val(k)./rad(k)) .* (xx(k)./rad(k)) .* (xx(k)./rad(k)) + val(k) ./ rad(k) );
            dv(k, 2) = coef *( (der(k) - val(k)./rad(k)) .* (xx(k)./rad(k)) .* (yy(k)./rad(k))                    );
            dv(k, 3) = coef *( (der(k) - val(k)./rad(k)) .* (xx(k)./rad(k)) .* (zz(k)./rad(k))                    );
            dv(~idxnz, 1) = coef * der(~idxnz);  % DX        
            countpp = countpp + 1;
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            vnlval(k) = coef *( (yy(k) ./ rad(k)) .* val(k) );
            dv(k, 1) = coef *( (der(k) - val(k)./rad(k)) .* (yy(k)./rad(k)) .* (xx(k)./rad(k))                    );
            dv(k, 2) = coef *( (der(k) - val(k)./rad(k)) .* (yy(k)./rad(k)) .* (yy(k)./rad(k)) + val(k) ./ rad(k) );
            dv(k, 3) = coef *( (der(k) - val(k)./rad(k)) .* (yy(k)./rad(k)) .* (zz(k)./rad(k))                    );
            dv(~idxnz, 2) = coef * der(~idxnz);  % DY        
            countpp = countpp + 1;
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            vnlval(k) = coef *( (zz(k) ./ rad(k)) .* val(k) );
            dv(k, 1) = coef *( (der(k) - val(k)./rad(k)) .* (zz(k)./rad(k)) .* (xx(k)./rad(k))                    );
            dv(k, 2) = coef *( (der(k) - val(k)./rad(k)) .* (zz(k)./rad(k)) .* (yy(k)./rad(k))                    );
            dv(k, 3) = coef *( (der(k) - val(k)./rad(k)) .* (zz(k)./rad(k)) .* (zz(k)./rad(k)) + val(k) ./ rad(k) );
            dv(~idxnz, 3) = coef * der(~idxnz);  % DZ
            countpp = countpp + 1;
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
        end
        
        if typ == PT.ptNLtype.L2
            % --------------------- d_z2 -------------------------
            coef = 1 / 4 * sqrt(5 * pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( - xx(k).*xx(k) - yy(k).*yy(k) + 2*zz(k).*zz(k) ) ./ (rad(k).*rad(k));
            Ylm(k, 2) = coef * ( -6 * xx(k) .* zz(k).^2 ./ rad(k).^4 );
            Ylm(k, 3) = coef * ( -6 * yy(k) .* zz(k).^2 ./ rad(k).^4 );
            Ylm(k, 4) = coef * (  6 * zz(k) .* (xx(k).^2 + yy(k).^2) ./ rad(k).^4 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;

            
            % -------------------- d_yz ----------------------------
            coef = 1 / 2 * sqrt(15 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( yy(k).*zz(k) ) ./ (rad(k).*rad(k));
            Ylm(k, 2) = coef * ( -2 * xx(k) .* yy(k) .* zz(k) ./ rad(k).^4 );
            Ylm(k, 3) = coef * (      zz(k) .* (zz(k).^2 + xx(k).^2 - yy(k).^2) ./ rad(k).^4 );
            Ylm(k, 4) = coef * (      yy(k) .* (yy(k).^2 + xx(k).^2 - zz(k).^2) ./ rad(k).^4 );
            
            vnlval(k, 1) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ------------------- d_xz -----------------------------
            coef = 1 / 2 * sqrt(15 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( zz(k).*xx(k) ) ./ (rad(k).*rad(k));
            Ylm(k, 2) = coef * (      zz(k) .* (zz(k).^2 - xx(k).^2 + yy(k).^2) ./ rad(k).^4 );
            Ylm(k, 3) = coef * ( -2 * xx(k) .* yy(k) .* zz(k) ./ rad(k).^4                   );
            Ylm(k, 4) = coef * (      xx(k) .* (xx(k).^2 + yy(k).^2 - zz(k).^2) ./ rad(k).^4 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ------------------- d_xy -----------------------------
            coef = 1 / 2 * sqrt(15 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( xx(k).*yy(k) ) ./ (rad(k).*rad(k));
            Ylm(k, 2) = coef * (      yy(k) .* (yy(k).^2 - xx(k).^2 + zz(k).^2) ./ rad(k).^4 );
            Ylm(k, 3) = coef * (      xx(k) .* (xx(k).^2 - yy(k).^2 + zz(k).^2) ./ rad(k).^4 );
            Ylm(k, 4) = coef * ( -2 * xx(k) .* yy(k) .* zz(k) ./ rad(k).^4                   );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;

            
            % ------------------- d_x^2-y^2 ------------------------
            coef = 1 / 4 * sqrt(15 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( xx(k).*xx(k) - yy(k).*yy(k) ) ./ (rad(k).*rad(k));
            Ylm(k, 2) = coef * (  2 * xx(k) .* (2*yy(k).^2 + zz(k).^2) ./ rad(k).^4 );
            Ylm(k, 3) = coef * ( -2 * yy(k) .* (2*xx(k).^2 + zz(k).^2) ./ rad(k).^4 );
            Ylm(k, 4) = coef * ( -2 * zz(k) .* (  xx(k).^2 - yy(k).^2) ./ rad(k).^4 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;

        end
        
        % FIXME: the derivative at r=0 for the f orbital MAY NOT BE CORRECT
        if typ == PT.ptNLtype.L3
            % ----------------------- f_z3 --------------------------------
            coef = 1 / 4 * sqrt(7 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( zz(k) .* (-3*xx(k).^2 - 3*yy(k).^2 + 2*zz(k).^2) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * (  3 * xx(k) .* zz(k)        .* (xx(k).^2 + yy(k).^2 - 4*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 3) = coef * (  3 * yy(k) .* zz(k)        .* (xx(k).^2 + yy(k).^2 - 4*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( -3 * (xx(k).^2 + yy(k).^2) .* (xx(k).^2 + yy(k).^2 - 4*zz(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ---------------------- f_y(3xx-yy) --------------------------
            coef = 1 / 4 * sqrt(35 / (2*pi));  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( yy(k) .* (3*xx(k).^2 - yy(k).^2) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * ( -3 * xx(k) .* yy(k) .* (xx(k).^2 - 3*yy(k).^2 - 2*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 3) = coef * (  3 * (xx(k).^4 - yy(k).^2 .* zz(k).^2 + xx(k).^2 .* (-3*yy(k).^2 + zz(k).^2)) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( -3 * yy(k) .* zz(k) .* (-3*xx(k).^2 + yy(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ---------------------- f_x(xx-3yy) --------------------------
            coef = 1 / 4 * sqrt(35 / (2*pi));  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( xx(k) .* (xx(k).^2 - 3*yy(k).^2) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * (  3 * (-yy(k).^2 .* (yy(k).^2 + zz(k).^2) + xx(k).^2 .* (3*yy(k).^2 + zz(k).^2)) ./ rad(k).^5 );
            Ylm(k, 3) = coef * (  3 * xx(k) .* yy(k) .* (-3*xx(k).^2 + yy(k).^2 - 2*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( -3 * xx(k) .* zz(k) .* (xx(k).^2 - 3*yy(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ---------------------- f_z(xx-yy) ---------------------------
            coef = 1 / 4 * sqrt(105 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( zz(k) .* (xx(k).^2 - yy(k).^2) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * ( xx(k) .* zz(k) .* (-xx(k).^2 + 5*yy(k).^2 + 2*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 3) = coef * ( yy(k) .* zz(K) .* (-5*xx(k).^2 + yy(k).^2 - 2*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( (xx(k).^2 - yy(k).^2) .* (xx(k).^2 + yy(k).^2 - 2*zz(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;

            % ----------------------- f_xyz -------------------------------
            coef = 1 / 2 * sqrt(105 / pi);  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( xx(k) .* yy(k) .* zz(k) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * ( yy(k) .* zz(k) .* (-2*xx(k).^2 + yy(k).^2 + zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 3) = coef * ( xx(k) .* zz(k) .* (xx(k).^2 - 2*yy(k).^2 + zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( xx(k) .* yy(k) .* (xx(k).^2 + yy(k).^2 - 2*zz(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ----------------------- f_yzz -------------------------------
            coef = 1 / 4 * sqrt(21 / (2*pi));  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( yy(k) .* (-xx(k).^2 - yy(k).^2 + 4*zz(k).^2) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * ( xx(k) .* yy(k) .* (xx(k).^2 + yy(k).^2 - 14*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 3) = coef * ( -(xx(k).^4 + 11*yy(k).^2 .*zz(k).^2 - 4*zz(k).^4 + xx(k).^2 .* (yy(k).^2 - 3*zz(k).^2)) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( yy(k) .* zz(k) .* (11*xx(k).^2 + 11*yy(k).^2 - 4*zz(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;
            
            % ----------------------- f_xzz -------------------------------
            coef = 1 / 4 * sqrt(21 / (2*pi));  % coefficients for spherical harmonics
            vnlval = zeros(idxsize, 1);  % value of nonlocal projectors
            dv = zeros(idxsize, dim);  % value of three derivatives
            Ylm = zeros(idxsize, dim+1);  % spherical harmonics (1) and its derivatives (2-4)
            idxnz = rad > MIN_RADIAL;
            k = idxnz;  % just use to simplify the code
            Ylm(k, 1) = coef * ( xx(k) .* (-xx(k).^2 - yy(k).^2 + 4*zz(k).^2) ) ./ (rad(k).*3);
            Ylm(k, 2) = coef * ( -(yy(k).^4 - 3*yy(k).^2 .*zz(k).^2 - 4*zz(k).^4 + xx(k).^2 .* (yy(k).^2 + 11*zz(k).^2)) ./ rad(k).^5 );
            Ylm(k, 3) = coef * ( xx(k) .* yy(k) .* (xx(k).^2 + yy(k).^2 - 14*zz(k).^2) ./ rad(k).^5 );
            Ylm(k, 4) = coef * ( xx(k) .* zz(k) .* (11*xx(k).^2 + 11*yy(k).^2 - 4*zz(k).^2) ./ rad(k).^5 );
            
            vnlval(k) = Ylm(k, 1) .* val(k);
            dv(k, 1) = Ylm(k, 1) .* der(k) .* (xx(k) ./ rad(k)) + Ylm(k, 2) .* val(k);
            dv(k, 2) = Ylm(k, 1) .* der(k) .* (yy(k) ./ rad(k)) + Ylm(k, 3) .* val(k);
            dv(k, 3) = Ylm(k, 1) .* der(k) .* (zz(k) ./ rad(k)) + Ylm(k, 4) .* val(k);
            
            countpp = countpp + 1; 
            vnlList.val(:, countpp) = vnlval;
            vnlList.drv{countpp} = dv;
            vnlList.wgt(countpp) = wgt;

        end    
    end
    
    % check the number of pseudpotentials
    if countpp ~= numpp
        error('countpp ~= numpp. Seriously wrong with nonlocal pseudopentials.');
    end
end

end
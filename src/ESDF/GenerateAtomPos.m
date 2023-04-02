function GenerateAtomPos(nreps, csize, atype, coefs, rands, file)
% GENERATEATOMPOS Generate the periodic atomic position via repeating given 
%    atom cell with potential randomness.
%
%    GenerateAtomPos(nreps, csize, natom, coefs, rands) writes the atom
%    positions after replication into file in reduced coordinate type.
%    input: 
%        nreps -- number of repeat in 3 direction, a 3-dim vector
%        csize -- size of single cell to perform periodic replication
%                 if csize is a scalar, treat all directions the same size 
%        atype -- list of types of each atom in a single cell
%        coefs -- reduced coordinates of atoms in a single cell
%        rands -- randomness add to atom position
%        file  -- file name to write atom position information
%                 default: atompos.txt
%
%    See also ESDF, ESDFReadInput.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


nrepx = nreps(1);
nrepy = nreps(2);
nrepz = nreps(3);
fprintf('Constructing periodic cells with repreat size %3d*%3d*%3d\n', ...
        nrepx, nrepy, nrepz);

natom = size(coefs, 1);
natomtot = natom * prod(nreps);
fprintf('Total number of atoms:   %6d\n\n', natomtot);

if numel(csize) == 1
    csize = ones(3, 1) * csize;
elseif numel(csize) == 3
    % do nothing
else
    error('csize is not in the correct form');
end
C = diag(csize);  % initial single cell

xyzmat = coefs * C';
atype = reshape(atype, [], 1);

%
% repeat the cell nrep times along the xyz directions
%
xyzlist = zeros(natomtot, 3);
typelist = zeros(natomtot, 1);
cnt = 0;
for krep = 1 : nrepz
    for jrep = 1 : nrepy
        for irep = 1 : nrepx
            xyzpick = xyzmat;  % atom position in single cell
            xyzpick(:, 1) = xyzpick(:, 1) + (irep-1) * csize(1);
            xyzpick(:, 2) = xyzpick(:, 2) + (jrep-1) * csize(2);
            xyzpick(:, 3) = xyzpick(:, 3) + (krep-1) * csize(3);
            idxsta = cnt * natom + 1;
            idxend = idxsta + natom - 1;
            xyzlist(idxsta : idxend, :) = xyzpick;
            typelist(idxsta : idxend) = atype;
            cnt = cnt + 1;
        end
    end
end

% Add randomness
xyzlist = xyzlist + rands * (randn(size(xyzlist)));

%
% modify the supercell
%
C(1, 1) = nrepx * csize(1);
C(2, 2) = nrepy * csize(2);
C(3, 3) = nrepz * csize(3);

fprintf('superCell C\n');
fprintf('%15.6f    %15.6f    %15.6f\n', C(1, 1), C(2, 2), C(3, 3));
fprintf('atom positions in (x, y, z) format\n');
fprintf('%15.6f    %15.6f    %15.6f\n', xyzlist');  % Transpose is important!

% 
% write atom positions to file
%
if nargin == 5
    fid = fopen('atompos.txt', 'w');
else
    fid = fopen(file, 'w');
end

tempC = diag(C);
tempC = reshape(tempC, 1, []);
fprintf(fid, 'begin Super_Cell\n');
fprintf(fid, '%12.6f    %12.6f    %12.6f\n', tempC(1), tempC(2), tempC(3));
fprintf(fid, 'end Super_Cell\n\n');

ntype = length(atype);
fprintf(fid, 'Atom_Types_Num:   %6d\n\n', ntype);

for i = 1 : length(atype)
    atomtype = atype(i);
    fprintf(fid, 'Atom_Type:        %6d\n\n', atomtype);

    xyzmat = xyzlist(typelist == atomtype, :);
    xyzmat = xyzmat ./ tempC;
    fprintf(fid, 'begin Atom_Red\n');
    fprintf(fid, '%12.6f    %12.6f    %12.6f\n', xyzmat');
    fprintf(fid, 'end Atom_Red\n\n');
end

fclose(fid);

end
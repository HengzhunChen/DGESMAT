%% test for LiNa

% generate atom position via replicating single cell
nreps = [1, 1, 4];
csize = 3.8507714854056379 / au2angDef();
atype = [3; 11];  % Li and Na
coefs = [
    0.0000000000000000    0.0000000000000000    0.0000000000000000
    0.5000000000000000    0.5000000000000000    0.5000000000000000
];
rands = 0;
GenerateAtomPos(nreps, csize, atype, coefs, rands);


% generate input file for pwdft
copyfile atompos.txt LiNa_pwdft.in

inputfile = 'LiNa_pwdft.in';
fid = fopen(inputfile, 'a');
fprintf(fid, 'SCF_Outer_MaxIter:    %d\n\n', 20);
fprintf(fid, 'SCF_Outer_Tolerance:  %12.6e\n\n', 1e-6);
fprintf(fid, 'Ecut_Wavefunction:    %12.6e\n\n', 20);
fprintf(fid, 'Pseudo_Type:          %s\n\n', 'ONCV');

outfile = 'statfile_pw';
pwdft_main(inputfile, outfile);


% generate input file for dgdft
copyfile atompos.txt LiNa_dgdft.in

inputfile = 'LiNa_dgdft.in';
fid = fopen(inputfile, 'a');
fprintf(fid, 'SCF_Outer_MaxIter:    %d\n\n', 20);
fprintf(fid, 'SCF_Outer_Tolerance:  %12.6e\n\n', 1e-6);
fprintf(fid, 'Ecut_Wavefunction:    %12.6e\n\n', 20);
fprintf(fid, 'Pseudo_Type:          %s\n\n', 'ONCV');
fprintf(fid, 'ALB_Num:              %d\n\n', 24);

fprintf(fid, 'begin Element_Size\n');
fprintf(fid, '%5d    %5d    %5d\n', 1, 1, 4);
fprintf(fid, 'end Element_Size\n\n');

outfile = 'statfile_dg';
dgdft_main(inputfile, outfile);

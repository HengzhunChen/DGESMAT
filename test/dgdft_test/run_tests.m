% test collection of DGDFT

% NOTE: Currently some tests use reference from package ScalES.
% Since cutoff value and RGaussian is different in two package,
% some tests may not be exactly the same.

count = 0;
failFlagList = [];

count = count + 1;
failFlagList(count) = test_default_value;

count = count + 1;
failFlagList(count) = test_mixing_density;

count = count + 1;
failFlagList(count) = test_mixing_kerker;

count = count + 1;
failFlagList(count) = test_periodic_potential;

count = count + 1;
failFlagList(count) = test_no_atom_density;

count = count + 1;
failFlagList(count) = test_hgh;

count = count + 1;
failFlagList(count) = test_vdw;

count = count + 1;
failFlagList(count) = test_xc_gga_pbe;

count = count + 1;
failFlagList(count) = test_buffer_size;


outputFid = fopen('test_report.out', 'w');
testdate = datetime;
fprintf(outputFid, 'Test date: %s \n', string(testdate));
fprintf(outputFid, '************************************************************\n');
for i = 1 : count
    if failFlagList(i) == 0
        fprintf(outputFid, 'test No. %d pass \n', i);
    else
        fprintf(outputFid, 'test No. %d fail \n', i);
    end
end
fprintf(outputFid, '************************************************************\n');
fclose(outputFid);

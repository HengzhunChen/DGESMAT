% test collection of DGDFT

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
failFlagList(count) = test_upf_atom_density;

count = count + 1;
failFlagList(count) = test_upf_oncv_pp;

count = count + 1;
failFlagList(count) = test_vdw;

count = count + 1;
failFlagList(count) = test_xc_gga_pbe;

count = count + 1;
failFlagList(count) = test_buffer_size;


outputFid = fopen('test_report.txt', 'w');
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

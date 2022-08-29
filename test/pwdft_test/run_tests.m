% test collection of PWDFT

count = 0;
failFlagList = [];

count = count + 1;
failFlagList(count) = test_default_value;

count = count + 1;
failFlagList(count) = test_dynamic_eigtol;

count = count + 1;
failFlagList(count) = test_mixing_kerker;

count = count + 1;
failFlagList(count) = test_pwsolver_chefsi;

count = count + 1;
failFlagList(count) = test_pwsolver_eigs;

count = count + 1;
failFlagList(count) = test_pwsolver_ppcg;

count = count + 1;
failFlagList(count) = test_upf_atom_density;

count = count + 1;
failFlagList(count) = test_upf_oncv_pp;

count = count + 1;
failFlagList(count) = test_use_vlocal;

count = count + 1;
failFlagList(count) = test_vdw;

count = count + 1;
failFlagList(count) = test_xc_gga_pbe;



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

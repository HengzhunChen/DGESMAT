function failFlag = test_upf_atom_density()
% test user option: isUseAtomDensity, i.e., use atom density from upf file
% as initial value of density

clear;
clear class;
rng(19981122);

timeStart = tic;

run('../../DGESMAT_startup.m');

tol = 1e-2;

checkListName = {...
    'H2O', ...
    'SiH4', ...
    };
checkList = zeros(length(checkListName), 1);
checkCount = 0;

%%
% Check for H2O molecule

inputFile = "./test_data/data_H2O/H2O_test_upf_rhoAtom.in";
outputFile = "./test_data/data_H2O/H2O_statfile";
info = pwdft_main(inputFile, outputFile);

Eref = -3.20446598e+00;

checkCount = checkCount + 1;
checkList(checkCount) = abs(info.Etot - Eref) / abs(Eref) < tol;
check_report(checkListName{checkCount}, Eref, info.Etot, checkList(checkCount));


%%
% Check for SiH4 molecule

inputFile = "./test_data/data_SiH4/SiH4_test_upf_rhoAtom.in";
outputFile = "./test_data/data_SiH4/SiH4_statfile";
info = pwdft_main(inputFile, outputFile);

Eref = -5.16435006e+00;

checkCount = checkCount + 1;
checkList(checkCount) = abs(info.Etot - Eref) / abs(Eref) < tol;
check_report(checkListName{checkCount}, Eref, info.Etot, checkList(checkCount));


%%
totTime = toc(timeStart);

fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('      Test for using UPF atom density to initialize \n');
fprintf('                     Overall Report\n');
fprintf('-------------------------------------------------------------\n');
for it = 1 : length(checkListName)
    fprintf('%18s   %d\n', checkListName{it}, checkList(it));
end
fprintf('-------------------------------------------------------------\n');
fprintf('%18s   %f seconds\n','Total Running Time:', totTime);
fprintf('=============================================================\n');

failFlag = sum(checkList == 0);

end

%%
function check_report(testName, Eref, E, check)
    fprintf('\n=============================================================\n');
    fprintf('Test Name:         %s\n', testName);
    fprintf('Reference Energy:  %s\n', Eref);
    fprintf('Current Energy:    %s\n', E);
    fprintf('Energy Accuracy:   %e\n', abs(Eref-E)/abs(Eref));
    if check
        fprintf('Test Result:      Passed\n');
    else
        fprintf('Test Result:      Failed\n');
    end
    fprintf('=============================================================\n');
end

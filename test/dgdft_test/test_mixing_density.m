function failFlag = test_mixing_density()
% test mixing with variable density

clear;
clear class;
rng(19981122);

timeStart = tic;

run('../../DGESMAT_startup.m');

tol = 1e-2;

checkListName = {...
    'H2', ...
    'LiNa', ...
    };
checkList = zeros(length(checkListName), 1);
checkCount = 0;

%%
% Check for H2 molecule (8 atoms)

inputFile = "./test_data/data_H2/H2_test_mixing_density.in";
outputFile = "./test_data/data_H2/H2_statfile";
info = dgdft_main(inputFile, outputFile);

Eref = -4.55239406e+00;

checkCount = checkCount + 1;
checkList(checkCount) = abs(info.Etot - Eref) / abs(Eref) < tol;
check_report(checkListName{checkCount}, Eref, info.Etot, checkList(checkCount));


%%
% Check for LiNa molecule (8 atoms)

inputFile = "./test_data/data_LiNa/LiNa_test_mixing_density.in";
outputFile = "./test_data/data_LiNa/LiNa_statfile";
info = dgdft_main(inputFile, outputFile);

Eref = -1.98239258e+02;

checkCount = checkCount + 1;
checkList(checkCount) = abs(info.Etot - Eref) / abs(Eref) < tol;
check_report(checkListName{checkCount}, Eref, info.Etot, checkList(checkCount));


%%
totTime = toc(timeStart);

fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('               Test for mixing with denisty\n');
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

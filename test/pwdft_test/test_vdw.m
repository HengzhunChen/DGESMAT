function failFlag = test_vdw()
% test with van der Waals energy and force.

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

inputFile = "./test_data/data_H2O/H2O_test_vdw.in";
outputFile = "./test_data/data_H2O/H2O_statfile";
info = pwdft_main(inputFile, outputFile);

Eref = -1.6809179687977204e-04;

checkCount = checkCount + 1;
checkList(checkCount) = abs(info.EVdw - Eref) / abs(Eref) < tol;
check_report(checkListName{checkCount}, Eref, info.EVdw, checkList(checkCount));


%%
% Check for SiH4 molecule

inputFile = "./test_data/data_SiH4/SiH4_test_vdw.in";
outputFile = "./test_data/data_SiH4/SiH4_statfile";
info = pwdft_main(inputFile, outputFile);

Eref = -2.4918862064051165e-03;

checkCount = checkCount + 1;
checkList(checkCount) = abs(info.EVdw - Eref) / abs(Eref) < tol;
check_report(checkListName{checkCount}, Eref, info.EVdw, checkList(checkCount));


%%
totTime = toc(timeStart);

fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('              Test for van der Waals energy\n');
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

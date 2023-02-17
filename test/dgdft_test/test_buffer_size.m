function failFlag = test_buffer_size()
% test with different buffer size, using reference from PWDFT with high 
% Ecut.
%
% buffer size is number of elements in the buffer region along +(-) x(y,z) 
% direction. Number of elements in extended element total will be 
%    (2 * bufferSize + 1) ^ d
% where d is number of dimensions that have been partitioned. 

clear;
clear class;
rng(19981122);

run('../../DGESMAT_startup.m');

timeStart = tic;

checkListName = {...
    'H2', ...
    };
checkList = zeros(length(checkListName), 1);
checkCount = 0;

%%
% Check for H2 molecule (8 atoms)

% compute energy reference using PWDFT
inputFile = "./test_data/data_H2/H2_PWDFT.in";
outputFile = "./test_data/data_H2/H2_statfile";
info = pwdft_main(inputFile, outputFile);
Eref = info.Etot; 

accuracyList = zeros(3, 1);

inputFile = "./test_data/data_H2/H2_test_buffer_size_1.in";
outputFile = "./test_data/data_H2/H2_statfile";
info = dgdft_main(inputFile, outputFile);
accuracyList(1) = abs(info.Etot - Eref) / abs(Eref);

inputFile = "./test_data/data_H2/H2_test_buffer_size_2.in";
outputFile = "./test_data/data_H2/H2_statfile";
info = dgdft_main(inputFile, outputFile);
accuracyList(2) = abs(info.Etot - Eref) / abs(Eref);

inputFile = "./test_data/data_H2/H2_test_buffer_size_3.in";
outputFile = "./test_data/data_H2/H2_statfile";
info = dgdft_main(inputFile, outputFile);
accuracyList(3) = abs(info.Etot - Eref) / abs(Eref);

checkCount = checkCount + 1;
if accuracyList(3) <= accuracyList(2) && accuracyList(2) <= accuracyList(1)
    checkList(checkCount) = 1;
else
    checkList(checkCount) = 0;
end
check_report(checkListName{checkCount}, Eref, accuracyList, checkList(checkCount));


%%
totTime = toc(timeStart);

fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('                Test for default values\n');
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
function check_report(testName, Eref, accuracyList, check)
    fprintf('\n=============================================================\n');
    fprintf('Test Name:         %s\n', testName);
    fprintf('Reference Energy:  %s\n', Eref);
    for i = 1 : length(accuracyList)
        fprintf('Accuracy for buffer size %d:  %e\n', i, accuracyList(i));
    end
    if check
        fprintf('Test Result:      Passed\n');
    else
        fprintf('Test Result:      Failed\n');
    end
    fprintf('=============================================================\n');
end

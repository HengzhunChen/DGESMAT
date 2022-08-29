function esdf_init(inputFile)
% ESDF_INIT reads which esdf parameters have been used in inputFile.
% 
%    Read keywords appeared in the inputFile into global variables lineList 
%    and tokenList. But not to initialize their value, concrete value is 
%    accessed by esdf_get(). This step is needed since labels in input file 
%    may be different from standard keyword list kw_label[].
%
%    See also ESDFReadInput, esdf_key, esdf_get.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global kw_label lineList tokenList

% define comment characters
comment = ['#', ';', '!'];
divide = {' ', '=', ':'};

% generate standard keyword list
esdf_unit();
esdf_key();

% "reduce" the keyword list for comparison
kw_label = esdf_reduce(kw_label);

% open the esdf input file
inputFid = fopen(inputFile, 'r');
% read the records from input file
i = 1;
while ~feof(inputFid)
    tline = fgetl(inputFid);
    % remove comment part
    for c = comment
        idxc = find(tline == c);
        if ~isempty(idxc)
            tline = tline(1: idxc(1)-1);
        end
    end
    % remove extra blank space
    if ~isempty(strtrim(tline))
        lineList(i) = strtrim(tline);
        i = i + 1;
    end
end
fclose(inputFid);

% read tokens, i.e. label and it value, from linelist
i = 1;
for line = lineList
    tokenLine = strsplit(line, divide, 'CollapseDelimiters', true);
    for j = 1 : length(tokenLine)
        tokenList(i, j) = tokenLine(j);
    end
    i = i + 1;
end

% data missing position will not use, here just for safety
tokenList(ismissing(tokenList)) = "";
% remove blank space
tokenList = erase(tokenList, ' ');
tokenList = erase(tokenList, '\t');


% check if any labels in the input file are unrecognized
temp_token = esdf_reduce(tokenList);
i = 1;
while i <= length(tokenList)
    % check if we are in a block
    label1 = temp_token(i, 1);
    if label1 == "begin"
        % check if block label is recognized
        label2 = temp_token(i, 2);
        idx_label2 = find(kw_label == label2, 1); 
        if isempty(idx_label2)
            error('Label "%s" not in keyword list', label2);
        end
        % check if block label is multiply defined
        if label2 ~= "atom_bohr" && label2 ~= "atom_ang" && label2 ~= "atom_red" 
            if sum(temp_token(:, 2) == label2) >= 4
                error('Label "%s" is multiply defined in the input file', label2);
            end
        end
        while temp_token(i, 1) ~= "end"
            i = i + 1;
        end
        i = i + 1;
    else  % not in block
        % check if label is in the list of keywords
        idx_label1 = find(kw_label == label1, 1);
        if isempty(idx_label1)
            error('Label "%s" not in keyword list', label1);
        end
        % check if label is multiply defined
        if label1 ~= "atom_type" && sum(temp_token(:, 1) == label1) >= 2
            error('Label "%s" is multiply defined in the input file', label1);
        end
        i = i + 1;
    end
end

% special check for atomList
numType_idx = temp_token(:, 1) == 'atom_types_num';
numAtomType = str2double(temp_token(numType_idx, 2));

if sum( temp_token(:, 1) == "atom_type" ) ~= numAtomType
    error('Atom type number is not match with atom data');
end
if sum( temp_token(:, 2) == "atom_bohr" ) ~= 0 && ...
   sum( temp_token(:, 2) == "atom_bohr" ) ~= 2*numAtomType
    error('Atom type number is not match with atom data');
elseif sum( temp_token(:, 2) == "atom_ang" ) ~= 0 && ...
       sum( temp_token(:, 2) == "atom_ang" ) ~= 2*numAtomType
    error('Atom type number is not match with atom data');
elseif sum( temp_token(:, 2) == "atom_red" ) ~= 0 && ...
       sum( temp_token(:, 2) == "atom_red" ) ~= 2*numAtomType
    error('Atom type number is not match with atom data');
end


end
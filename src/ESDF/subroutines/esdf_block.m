function [block_data, nlines] = esdf_block(block_label, atom_type_index)
% ESDF_BLOCK reads block data and number of lines attached to block label.
%
%    [block_data, nlines] = esdf_block(block_label) returns block data and
%    number of lines corresponding to the block_label, block_data is a 
%    column vector in types of string. If block_data is not found, returns
%    "" and 0.
%
%    [block_data, nlines] = esdf_block(block_label, atom_type_index)
%    atom_type_index is used only when read block of atom position. 
%
%    See also ESDFReadInput. 

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global kw_label tokenList lineList;

% check label is defined
labelstr = esdf_reduce(block_label);
index = find(kw_label == labelstr, 1);
if isempty(index)
    error('Label "%s" not recognized in keyword list', block_label);
end

temp_token = esdf_reduce(tokenList);

if nargin == 1
    % read data from input file by tokenList
    block_data = "";
    nlines = 0;
    token_idx = find(temp_token(:, 1) == "begin");
    for idx = token_idx'
        if temp_token(idx, 2) == esdf_reduce(block_label)
            i = 1;
            while temp_token(idx+i, 1) ~= "end"
                % note: we use lineList here
                block_data(i) = lineList(idx + i);
                block_data(i) = replace(block_data(i), ',', ' ');
                block_data(i) = replace(block_data(i), ';', ' ');
                i = i + 1;
            end
            nlines = i - 1;
        end
    end
elseif nargin == 2
    % read data of atom position
    block_data = "";
    nlines = 0;
    count_atom = 1;
    token_idx = find(temp_token(:, 1) == "begin");
    for idx = token_idx'
        if temp_token(idx, 2) == esdf_reduce(block_label)
            if count_atom == atom_type_index
                i = 1;
                while temp_token(idx+i, 1) ~= "end"
                    block_data(i) = lineList(idx + i);
                    block_data(i) = replace(block_data(i), ',', ' ');
                    block_data(i) = replace(block_data(i), ';', ' ');
                    i = i + 1;
                end
                nlines = i - 1;
                break;
            else
                count_atom = count_atom + 1;
                continue;
            end
        end
    end
end
 
end
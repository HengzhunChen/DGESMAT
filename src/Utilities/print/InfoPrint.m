function InfoPrint(varargin)
% INFOPRINT prints variables and their value. 
% 
%    InfoPrint(msg) prints message to statistics file, screen and debug 
%    file.
%    
%    InfoPrint(name, value) prints variable with name and value to 
%    statistics file, screen and debug file. 
%
%    InfoPrint(name, value, unit) prints variable and its value and unit to
%    statistics file, screen and debug file.
%
%    InfoPrint(nameList, valueList) prints a list to variables and their
%    values to statistics file, screen and debug file.
%
%    InfoPrint(formatSpec, A1, ..., An) applies the formatSpec to all 
%    elements of arrays A1,...An in column order, and writes the data to 
%    statistics file, screen and debug file. 
%
%
%    InfoPrint(printId, msg) prints message to output files with ID 
%    printId.
%    
%    InfoPrint(printId, name, value) prints variable with name and value to 
%    output file with ID printId. 
%
%    InfoPrint(printId, name, value, unit) prints variable and its value 
%    and unit to output file with ID printId.
%
%    InfoPrint(printId, nameList, valueList) prints a list to variables 
%    and their values to output file with ID printId.
%
%    InfoPrint(printId, formatSpec, A1, ..., An) applies the formatSpec to 
%    all elements of arrays A1,...An in column order, and writes the data 
%    to output file with ID printId. 
%
%    If printId == 0, print the banner to the statistics file.
%    If printId == 1, print the banner to the screen.
%    If printId == 2, print the banner to debug file.
%
%    See also InfoPrint, PrintBlock, fprintf.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global outFid debugFid;
screenID = 1;

printId = []; 
if isa(varargin{1}, 'numeric')
    fidList = varargin{1};
    infoStart = 2;
    if ismember(0, fidList)
        printId = [printId, outFid];
    end
    if ismember(1, fidList)
        printId = [printId, screenID];
    end
    if ismember(2, fidList)
        printId = [printId, debugFid];
    end
else
    infoStart = 1;
    printId = [outFid, screenID, debugFid];
end

name = varargin{infoStart};
if contains(name, '%') || contains(name, '\')
    if nargin == infoStart
        for fid = printId
            fprintf(fid, name);
        end
    else
        valueList = varargin(infoStart+1 : end);
        if length(valueList) == 1
            value = valueList{1};
        else
            value = cell2mat(valueList);
        end
        for fid = printId
            fprintf(fid, name, value);
        end
    end
else
    if nargin == infoStart
        name = varargin{infoStart};
        for fid = printId
            fprintf(fid, '%s\n', name);
        end
    elseif nargin == infoStart + 1
        name = varargin{infoStart};
        value = varargin{infoStart + 1};
        if isnumeric(value)
            if contains(num2str(value), '.')
                valuestr = num2str(value, '%+1.8e ');
            else
                valuestr = num2str(value);
            end
        else
            valuestr = string(value);
        end
        for fid = printId
            fprintf(fid, '%s%s\n', name, valuestr);
        end
    elseif nargin == infoStart + 2
        name = varargin{infoStart};
        value = varargin{infoStart + 1};
        unit = varargin{infoStart + 2};

        if isnumeric(value)
            if contains(num2str(value), '.')
                valuestr = num2str(value, '%+1.8e ');
            else
                valuestr = num2str(value);
            end
        else 
            valuestr = string(value);
        end
        for fid = printId
            fprintf(fid, '%s%s %s\n', name, valuestr, unit);
        end
    else
        nameList = string( varargin( (infoStart) : 2 : end ) );
        valueList = varargin( (infoStart+1) : 2 : end );
        valuestr = cell(size(valueList));

        for i = 1 : length(valueList)
            value = valueList{i};
            if isnumeric(value)
                if contains(num2str(value), '.')
                    valuestr{i} = num2str(value, '%+1.8e ');
                else
                    valuestr{i} = num2str(value);
                end
            else 
                valuestr{i} = string(value);
            end
        end

        numVar = nargin - 1;
        for fid = printId
            for i = 1 : numVar/2
                fprintf(fid, '%s %s   ', nameList(i), valuestr{i});
            end
            fprintf(fid, '\n');
        end
    end
end

end
function PrintBlock(varargin)
% PRINTBLOCK print the message in a block banner.
%
%    PrintBlock(msg) prints the message in the banner to statistics file
%    and screen.
%
%    PrintBlock(variableName, value) prints the variable with variableName
%    and value in the banner to statistics file and screen.
%
%    PrintBlock(printId, msg) prints the message in the banner to output 
%    files with ID printId.
% 
%    PrintBlock(printId, variableName, value) prints the variable with 
%    value in the banner to output files with ID printId.
%
%    If printId == 0, print the banner to the statistics file.
%    If printId == 1, print the banner to the screen.
%
%    See also InfoPrint, fprintf.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


global outFid;
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
else
    infoStart = 1;
    printId = [outFid, screenID];
end

if nargin == infoStart
    name = varargin{infoStart};
    for fid = printId
        fprintf(fid, '\n*********************************************************************\n');
        fprintf(fid, '%s \n', name);
        fprintf(fid, '*********************************************************************\n\n');
    end
elseif nargin == infoStart + 1
    name = varargin{infoStart};
    value = varargin{infoStart + 1};
    valuestr = string(value);
    for fid = printId
        fprintf(fid, '\n*********************************************************************\n');
        fprintf(fid, '%s%s \n', name, valuestr);
        fprintf(fid, '*********************************************************************\n\n');
    end
end    

end
function define_global(outFile)
% DEFINE_GLOBAL define global variables
%
%    define_global(outFile) initializes some global variables 
%    and is used as a reference for the usage of the global variables.
%
%    See also dgdft_main, pwdft_main.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


%
% output files
%

global outFid
outFid = fopen(outFile, 'w');
% used in
% pwdft_main()
% dgdft_main()
% InfoPrint()
% PrintBlock()


%
% global variables used in esdfParam to read input
%

% global variables deal with conversion of different physical units
global phy_var phy_name phy_unit;
phy_var = "";
phy_name = "";
phy_unit = [];
% used in esdf_unit(), esdf_convfac()

global kw_label lineList tokenList;
% kw_label is the list of all ESDF parameters
% lineList is the line list read from input file, used for extracting
% tokenList and read data for block type.
% tokenList is the list of tokens, each row is label and its value, not
% including the block type labels.
kw_label = "";
lineList = "";
tokenList = "";
% used in esdf_init(), esdf_key(), esdf_block(), esdf_get()


end
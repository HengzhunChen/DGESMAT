function varargout = dgdft_main(varargin)
% DGDFT_MAIN main function of dgdft
%    Read parameters of dgdft from inputFile, then run the routines of
%    DGDFT and print outcome in outFile.
%
%    dgdft(inputFile) runs the routines of DGDFT according to the data from
%    file inputFile and print the outcome to default output file statfile.
%
%    dgdft(inputFile, outFile) runs the routines of DGDFT according to the
%    data from file inputFile and print the outcome to file outFile.
%
%    dgdft(inputFile, outFile, debugFile) runs the routines of DGDFT
%    according to data from file inputFile and print the outcome to file
%    outFile, debugFile is used when you set the print ID in InfoPrint() 
%    equals to 2 and prints some data you need into file debugFile.
%
%    info = dgdft(inputFile, __, __) runs the routines of DGDFT and output
%    some information into struct info, which contains variables Etot, 
%    Evdw and Efree. 
%
%    See also InfoPrint.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


if nargin == 1
    inputFile = varargin{1};
    outFile   = "statfile";
    debugFile = "";
elseif nargin == 2
    inputFile = varargin{1};    
    outFile   = varargin{2};
    debugFile = "";
elseif nargin == 3
    inputFile = varargin{1};
    outFile   = varargin{2};
    debugFile = varargin{3};
else
    error('Wrong number of arguments');
end

% setup global variables
define_global(outFile, debugFile);
global outFid debugFid esdfParam;

InfoPrint(1, "run dgdft \n");
InfoPrint(1, "inFile:    %s \n", inputFile);


% =======================================================================
% Input parameters
% =======================================================================

% Initialize input parameters
timeStart = tic;
ESDFReadInput( inputFile ); 
timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for reading the input file is %8f [s]\n\n', timeEnd);

ESDFPrintInput();


% =======================================================================
% Preparation
% =======================================================================

% Setup periodic table
timeStart = tic;
ptable = PeriodTable();
timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for setting up the periodic table is %8f [s]\n\n', timeEnd);

% Setup the element and extended element information
VecEigSol = setup_element(esdfParam, ptable);

% Setup HamiltonianDG
timeStart = tic;
hamDG = HamiltonianDG();
hamDG = CalculatePseudoPotential(hamDG, ptable);
timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for setting up the DG Hamlitonian is %8f [s]\n\n', timeEnd);

% Setup SCFDG
timeStart = tic;
scfDG = SCFDG(hamDG, VecEigSol, ptable);
timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for setting up SCFDG is %8f [s]\n\n', timeEnd);


% ======================================================================
% Solve: a single shot calculation
% ======================================================================

% Main SCF iteration
timeStart = tic;
scfDG = Iterate(scfDG);
timeEnd = toc(timeStart);
InfoPrint([0, 1], '! Total time for the SCF iteration = %8f [s]\n\n', timeEnd);

InfoPrint(1, 'End of single dgdft\n');
InfoPrint(1, 'outFile:    %s \n', outFile);

fclose(outFid);
if ~isempty(debugFid)
    fclose(debugFid);
end


% ======================================================================
% TODO: Geometry optimization or Molecular dynamics
% ======================================================================



% ======================================================================
% output part of information if needed
% ======================================================================

if nargout == 1
    info.Etot = scfDG.Etot;
    info.Evdw = scfDG.Evdw;
    info.Efree = scfDG.Efree;
    
    varargout{1} = info;
end


end
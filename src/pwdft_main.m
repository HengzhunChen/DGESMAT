function varargout = pwdft_main(varargin)
% PWDFT_MAIN main function of pwdft
%    Read parameters of PWDFT from inputFile, then run the routines of 
%    PWDFT and print data and outcome to outFile.
%
%    pwdft(inputFile) runs the routines of PWDFT according to the data from
%    file inputFile and print the outcome to default output file statfile.
%
%    pwdft(inputFile, outFile) runs the routines of PWDFT according to the
%    data from file inputFile and print the outcome to file outFile.
%
%    pwdft(inputFile, outFile, debugFile) runs the routines of PWDFT
%    according to data from file inputFile and print the outcome to file
%    outFile, debugFile is used when you set the print ID in InfoPrint() 
%    equals to 2 and prints some data you need into file debugFile.
%
%    info = pwdft(inputFile, __, __) runs the routines of PWDFT and output
%    some information into struct info, which contains variables Etot, 
%    Evdw and Efree. 
%
%    See also InfoPrint.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


if nargin == 1
    inputFile = varargin{1};
    outFile = "statfile";
    debugFile = "";
elseif nargin == 2
    inputFile = varargin{1};
    outFile = varargin{2};
    debugFile = "";
elseif nargin == 3
    inputFile = varargin{1};
    outFile = varargin{2};
    debugFile = varargin{3};
else
    error('Wrong number of arguments');
end

% setup global variables
define_global(outFile, debugFile);
global outFid debugFid esdfParam;

InfoPrint(1, "run pwdft \n");
InfoPrint(1, "inFile:    %s \n", inputFile);


%========================================================================
%                          Input parameter
%========================================================================

% Initialize input parameters
timeStart = tic;
ESDFReadInput( inputFile ); 
timeEnd = toc( timeStart );
InfoPrint([0, 1], 'Time for reading the input file is %8f [s]\n', timeEnd);

ESDFPrintInput();


%========================================================================
%                            Preparation
%========================================================================

% Hamiltonian
hamKS = HamiltonianKS(esdfParam);
timeStart = tic;
hamKS = CalculatePseudoPotential(hamKS);
timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for calculating the pseudopotential for the Hamlitonian is %8f [s]\n', timeEnd);

% Wavefunctions
domain = esdfParam.basic.domain;
psi = Spinor(domain, hamKS.NumStateTotal(), 0.0);
psi.wavefun = rand( size(psi.wavefun) );
% used for debug
% for k = 1 : hamKS.NumStateTotal()
%    psi.wavefun(k, k) = 1;
% end

% EigenSolverKS
eigSol = EigenSolverKS(hamKS, psi);

% SCF
scf = SCF(eigSol);


%=======================================================================
%                   Solve: a single shot calculation
%=======================================================================

% Main SCF iteration
timeStart = tic;
scf = Iterate(scf);
timeEnd = toc(timeStart);
InfoPrint([0, 1], '! Total time for the SCF iteration = %8f [s]\n', timeEnd);

InfoPrint(1, 'end of single pwdft \n');
InfoPrint(1, 'outFile:    %s \n', outFile);

fclose(outFid);
if ~isempty(debugFid)
    fclose(debugFid);
end


%=======================================================================
% Geometry optimization or Molecular dynamics
%=======================================================================

% TODO


% ======================================================================
% output part of information if needed
% ======================================================================

if nargout == 1
    info.Etot = scf.Etot;
    info.EVdw = scf.EVdw;
    info.Efree = scf.Efree;
    
    varargout{1} = info;
end


end
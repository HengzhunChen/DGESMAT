function PT = Setup(PT, esdfParam)
% PERIODTABLE/SETUP prepares data of PeriodTable class.
%
%    See also PeriodTable, PTEntry, PeriodChart, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


pseudoType = esdfParam.basic.pseudoType;
upfFile = esdfParam.basic.upfFile;
atomTypes = esdfParam.basic.atomTypes;

periodChart = PeriodChart();  % basic information of atoms

% read pseudopotential data
if pseudoType == "ONCV"
    if ~isempty(upfFile)
        % read data from UPF files
        for i = 1 : length(upfFile)
            [tempEntry, Znuc] = ReadUPF(upfFile{i});
            PT.pteMap(Znuc) = tempEntry;
        end
    else
        % use data from default UPF files
        for i = 1 : length(atomTypes)
            atomType = atomTypes(i);
            atomSymbol = periodChart.symbol(atomType);
            
            fileName = atomSymbol + "_ONCV_PBE-1.2.upf";
            upfFile = dgesmat_root() + ...
                "/ppdata/sg15_oncv_upf_2020/" + fileName;
            if ~exist(upfFile, 'file')
                fileName = atomSymbol + "_ONCV_PBE-1.1.upf";
                upfFile = dgesmat_root() + ...
                    "/ppdata/sg15_oncv_upf_2020/" + fileName;                
                if ~exist(upfFile, 'file')
                    fileName = atomSymbol + "_ONCV_PBE-1.0.upf";
                    upfFile = dgesmat_root() + ...
                     "/ppdata/sg15_oncv_upf_2020/" + fileName;
                end
            end

            [tempEntry, Znuc] = ReadUPF(upfFile);
            PT.pteMap(Znuc) = tempEntry;            
        end
    end
else 
    % generate pseudopotential data instead of reading from UPF file
    InfoPrint([0, 1], 'Generate HGH type pseudopotential as default. \n');

    % TODO: currently not support
    if PT.userOption.isUseVLocal
        error('HGH does not support option isUseVLocal currently');
    end
    Znucs = atomTypes;
    for i = 1 : length(Znucs)
        Znuc = Znucs(i);
        tempEntry = HGH(Znuc);
        PT.pteMap(Znuc) = tempEntry;
    end
end


% Extra processing of Vlocal data
if PT.userOption.isUseVLocal
    for key = PT.pteMap.keys
        type = key{1};

        Zion = PT.pteMap(type).params.Zion;
        RGaussian = PT.pteMap(type).params.RGaussian;
        rad = PT.pteMap(type).samples.rad;
        vlocal = PT.pteMap(type).samples.vLocal;
        
        % remove the pseudocharge contribution
        idxz = rad == 0;
        vlocal(idxz) = vlocal(idxz) + Zion / RGaussian * 2 / sqrt(pi);
        idxnz = rad ~= 0;
        vlocal(idxnz) = vlocal(idxnz) + ...
            ( Zion ./ rad(idxnz) ) .* erf(rad(idxnz) / RGaussian);
        
        ptEntry = PT.pteMap(type);
        ptEntry.samples.vLocal = vlocal;
        PT.pteMap(type) = ptEntry;
    end
end


% create splines
for key = PT.pteMap.keys
    type = key{1};
    
    splineTemp = struct(...
        'pseudoCharge', [], ...
        'drv_pseudoCharge', [], ...
        'vLocal', [], ...
        'drv_vLocal', [], ...
        'rhoAtom', [], ...
        'drv_rhoAtom', [], ...
        'nonlocal', [] ...
        );
    rad = PT.pteMap(type).samples.rad;
    
    if ~PT.userOption.isUseVLocal
        a = PT.pteMap(type).samples.pseudoCharge;
        [b, c, d] = spline(rad, a);
        splineTemp.pseudoCharge = [rad, a, b, c, d];

        a = PT.pteMap(type).samples.drv_pseudoCharge;
        [b, c, d] = spline(rad, a);
        splineTemp.drv_pseudoCharge = [rad, a, b, c, d];
    else    
        a = PT.pteMap(type).samples.vLocal;
        [b, c, d] = spline(rad, a);
        splineTemp.vLocal = [rad, a, b, c, d];

        a = PT.pteMap(type).samples.drv_vLocal;
        [b, c, d] = spline(rad, a);
        splineTemp.drv_vLocal = [rad, a, b, c, d];
    end
    
    if esdfParam.userOption.general.isUseAtomDensity
        a = PT.pteMap(type).samples.rhoAtom;
        [b, c, d] = spline(rad, a);
        splineTemp.rhoAtom = [rad, a, b, c, d];

        a = PT.pteMap(type).samples.drv_rhoAtom;
        [b, c, d] = spline(rad, a);
        splineTemp.drv_rhoAtom = [rad, a, b, c, d];
    end

    NLsamples = PT.pteMap(type).samples.nonlocal;
    NLsize = size(NLsamples);
    splineTemp.nonlocal = zeros(NLsize(1), 5*NLsize(2));
    for i = 1 : NLsize(2)
        a = NLsamples(:, i);
        [b, c, d] = spline(rad, a);
        idxStart = 5*(i-1) + 1;
        idxEnd = idxStart + 4;
        splineTemp.nonlocal(:, idxStart : idxEnd) = [rad, a, b, c, d];
    end

    PT.splineMap(type) = splineTemp;
end


end
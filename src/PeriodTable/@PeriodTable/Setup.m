function PT = Setup(PT, esdfParam)
% PERIODTABLE/SETUP prepares data of PeriodTable class.
%
%    See also PeriodTable, PTEntry, ElementTable, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


pseudoType = esdfParam.basic.pseudoType;
upfFile = esdfParam.basic.upfFile;
atomTypes = esdfParam.basic.atomTypes;
domain = esdfParam.basic.domain; 

etable = ElementTable();  % basic information of atoms

%
% read pseudopotential data
%
if ~isempty(upfFile)
% if UPF file is provided in input file, use UPF file from input file
    for i = 1 : length(upfFile)
        [tempEntry, Znuc] = ReadUPF(upfFile{i});
        PT.pteMap(Znuc) = tempEntry;
    end
else
% if UPF file is not provided while pseudoType is assigned in input file,
% use default UPF file corresponding to pseudoType   
    if pseudoType == "ONCV"
        % use data from default UPF files
        InfoPrint([0, 1], 'Read default ONCV type pseudopotential. \n');
        ONCVpath = dgesmat_root() + "/ppdata/sg15_oncv_upf_2020/";
        for i = 1 : length(atomTypes)
            atomType = atomTypes(i);
            atomSymbol = etable.symbol(atomType);
            
            fileName = atomSymbol + "_ONCV_PBE-1.2.upf";
            upfFile = ONCVpath + fileName;
            if ~exist(upfFile, 'file')
                fileName = atomSymbol + "_ONCV_PBE-1.1.upf";
                upfFile = ONCVpath + fileName;                
                if ~exist(upfFile, 'file')
                    fileName = atomSymbol + "_ONCV_PBE-1.0.upf";
                    upfFile = ONCVpath + fileName;
                end
            end

            [tempEntry, Znuc] = ReadUPF(upfFile);
            PT.pteMap(Znuc) = tempEntry;            
        end
    elseif pseudoType == "HGH"
        InfoPrint([0, 1], 'Read default HGH type pseudopotential. \n');
        HGHpath = dgesmat_root() + "/ppdata/hgh_pbe_upf/";
        for i = 1 : length(atomTypes)
            atomType = atomTypes(i);
            atomSymbol = etable.symbol(atomType);
            
            upfFile = HGHpath + atomSymbol + ".pbe-hgh.UPF";
            if ~exist(upfFile, 'file')
                files = dir(HGHpath);
                names = {files.name};
                pattern = atomSymbol + ".";
                foundFlag = false;
                for j = 1 : length(names)
                    fileName = names{j};
                    if contains(fileName, pattern)
                        upfFile = HGHpath + fileName;
                        foundFlag = true;
                        break;
                    end
                end
                if ~foundFlag
                    msg = "There is no default UPF file for chemical element " ...
                        + atomSymbol + ". Please try other type pseudopotential " + ...
                        + "or provide your UPF file in input file.";
                    error(msg);
                end
            end

            [tempEntry, Znuc] = ReadUPF(upfFile);
            PT.pteMap(Znuc) = tempEntry;            
        end
    else
        error('pseudopotential type %s is not supported\n', pseudoType);
    end
end

% check for whether cutoffs greater than supercell size
dmlength = domain.length;
for i = 1 : length(atomTypes)
    atomType = atomTypes(i);
    vsrcut = RcutVLocalSR(PT, atomType);
    exceedSR = find(dmlength < 2*vsrcut);
    if ~isempty(exceedSR)
        msg = "short range cutoff of VLocal of chemical element " + etable.symbol(atomType) + ...
            " is greater than size of supercell in " + num2str(length(exceedSR)) + " directions." + ...
            " This will cause extra cutoff in current implementation.";
        warning(msg);
    end
    vnlcut = RcutNonlocal(PT, atomType);
    exceedNL = find(dmlength < 2*vnlcut);
    if ~isempty(exceedNL)
        msg = "nonlocal potential cutoff of chemical element " + etable.symbol(atomType) + ...
            " is greater than size of supercell in " + num2str(length(exceedNL)) + " directions." + ...
            " This will cause extra cutoff in current implementation.";
        warning(msg);
    end
    if esdfParam.userOption.general.isUseAtomDensity
        rhoatomcut = RcutRhoAtom(PT, atomType);
        exceedRA = find(dmlength < 2*rhoatomcut);
        if ~isempty(exceedRA)
            msg = "cutoff of rhoAtom (used as initial density) of chemical element " + etable.symbol(atomType) + ...
                " is greater than size of supercell in " + num2str(length(exceedRA)) + " directions." + ...
                " This will cause extra cutoff in current implementation.";
            warning(msg);
        end
    end
end


% Post-processing of VLocal data
% Seperate the VLocal into short range part VLocalSR and long range part VGaussian
for key = PT.pteMap.keys
    type = key{1};
    Zval = PT.pteMap(type).params.Zval;
    RGaussian = PT.pteMap(type).params.RGaussian;
    rad = PT.pteMap(type).samples.rad;
    vlocal = PT.pteMap(type).samples.vLocalSR;
    
    % remove the pseudocharge contribution
    vlocalSR = vlocal + VGaussian(rad, Zval, RGaussian);
    
    ptEntry = PT.pteMap(type);
    ptEntry.samples.vLocalSR = vlocalSR;
    PT.pteMap(type) = ptEntry;
end


% create splines
for key = PT.pteMap.keys
    type = key{1};
    
    splineTemp = struct(...
        'vLocalSR', [], ...
        'drv_vLocalSR', [], ...
        'rhoAtom', [], ...
        'drv_rhoAtom', [], ...
        'nonlocal', [] ...
        );
    rad = PT.pteMap(type).samples.rad;
    
    a = PT.pteMap(type).samples.vLocalSR;
    [b, c, d] = spline(rad, a);
    splineTemp.vLocalSR = [rad, a, b, c, d];

    a = PT.pteMap(type).samples.drv_vLocalSR;
    [b, c, d] = spline(rad, a);
    splineTemp.drv_vLocalSR = [rad, a, b, c, d];
    
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


function vGaussian = VGaussian(rad, Zval, RGaussian)
% VGAUSSIAN returns the potential due to a Gaussian compensation charge.    
    EPS = 1e-12;
    idxz = rad < EPS;
    idxnz = ~idxz;
    vGaussian = zeros(size(rad));
    vGaussian(idxnz) = Zval ./ rad(idxnz) .* erf(rad(idxnz) ./ RGaussian);
    vGaussian(idxz) = Zval ./ RGaussian * 2.0 ./ sqrt(pi);
end

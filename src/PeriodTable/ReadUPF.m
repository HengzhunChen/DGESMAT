function [tempEntry, atomicNum] = ReadUPF(fileName)
% READUPF reads pseudopotential info from UPF file
%
%    [tempEntry, atomicNum] = ReadUPF(fileName) reads data from file and
%    returns the atomic number, pseudopotential data is stored in the 
%    structure PTEntry.
%
%    See also PeriodTable, PTEntry.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

% NOTE:
% Currently only version 2 UPF file with norm-conserving condition can be
% supported.

pchart = PeriodChart();
tempEntry = PTEntry();

upfid = fopen(fileName, 'r');

% determine UPF version
upf_version = 0;

% The first line of the UPF potential file contains either of the
% following:
%    <PP_INFO> (for UPF version 1)
%    <UPF version="2.0.1"> (for UPF version 2)
% Here only version 2 is supported

buff = fgetl(upfid);
if contains(buff, "<PP_INFO>")
    upf_version = 1;
elseif contains(buff, '<UPF version="2.0.1">')
    upf_version = 2;
end
if upf_version == 0
    InfoPrint(0, ' Format of UPF file not recognized \n');
    InfoPrint(0, ' First line of file: %s \n', buff);
    error(' Format of UPF file not recognized ');
end
InfoPrint(0, 'UPF version: %d \n', upf_version);


if upf_version == 1
    error(' Format of UPF file 1.0 not supported');
elseif upf_version == 2
    % process UPF version 2 potential 
    
    % ---------------------- <PP_INFO> --------------------------------
    % read lines until <PP_INFO> is found
    while ~contains(buff, '<PP_INFO>')
        buff = fgetl(upfid);
    end
    upf_pp_info = "";
    while ~contains(buff, '</PP_INFO>')
        buff = fgets(upfid);
        upf_pp_info = upf_pp_info + buff; 
    end
    
    % remove all '<' and '>' characters from the PP_INFO field for XML
    % compatibility
    upf_pp_info = erase(upf_pp_info, '<');
    upf_pp_info = erase(upf_pp_info, '>');
    
    
    % -------------------- <PP_HEADER> ------------------------------
    tag = read_start_field("PP_HEADER", upfid);
    
    % get attribure "element"
    upf_symbol = get_attr(tag, "element");
    upf_symbol = erase(upf_symbol, ' ');
    
    % get atomic number and mass
    atomicNum = pchart.z(upf_symbol);
    mass = pchart.mass(upf_symbol);
    
    tempEntry.params.Znuc = atomicNum;
    tempEntry.params.Mass = mass;
    
    % check if potential is norm-conserving or semi-local
    pseudo_type = get_attr(tag, "pseudo_type");
    InfoPrint(0, ' pseudo_type = %s \n', pseudo_type);
    if pseudo_type ~= "NC" && pseudo_type ~= "SL"
        error('pseudo_type must be NC or SL');
    end
    
    % NLCC flag
    upf_nlcc_flag = get_attr(tag, "core_correction");
    if upf_nlcc_flag == "T"
        InfoPrint(0, ' Potential includes a non-linear core correction \n');
    end
    
    % XC functional 
    upf_functional = get_attr(tag, "functional");
    InfoPrint(0, ' upf_functional = %s \n', upf_functional);
    
    % valence charge
    upf_zval = get_attr(tag, "z_valence");
    InfoPrint(0, ' upf_zval = %s \n', upf_zval);
    tempEntry.params.Zion = str2double(upf_zval);
        
    % max angular momentum
    upf_lmax = get_attr(tag, "l_max");
    InfoPrint(0, ' upf_lmax = %s \n', upf_lmax);
    upf_lmax = str2double(upf_lmax);
    
    % local angualr moment
    upf_llocal = get_attr(tag, "l_local");
    InfoPrint(0, ' upf_llocal = %s \n', upf_llocal);
    upf_llocal = str2double(upf_llocal);
    
    % number of points in mesh
    upf_mesh_size = get_attr(tag, "mesh_size");
    InfoPrint(0, ' upf_mesh_size = %s \n', upf_mesh_size);
    upf_mesh_size = str2double(upf_mesh_size);
    
    % number of wavefunctions
    upf_nwf = get_attr(tag, "number_of_wfc");
    InfoPrint(0, ' upf_nwf = %s \n', upf_nwf);
    upf_nwf = str2double(upf_nwf);
    
    % number of projectors
    upf_nproj = get_attr(tag, "number_of_proj");
    InfoPrint(0, ' upf_proj = %s \n', upf_nproj);
    upf_nproj = str2double(upf_nproj);


    % FIXME: RGaussian should be given by a table according to the element type 
    tempEntry.params.RGaussian = 1.0;    
    tempEntry.params.Eself = 1.0;    
    
    tempEntry.samples.nonlocal = zeros(upf_mesh_size, 2*upf_nproj);
    tempEntry.cutoffs.nonlocal = zeros(upf_nproj, 1);
    tempEntry.NLweights        = zeros(upf_nproj, 1);
    tempEntry.NLtypes          = zeros(upf_nproj, 1);
    
    tempEntry.samples.drv_pseudoCharge = zeros(upf_mesh_size, 1);
    tempEntry.samples.drv_vLocal       = zeros(upf_mesh_size, 1);
    tempEntry.samples.drv_rhoAtom      = zeros(upf_mesh_size, 1);
    
    
    % TODO: the following rhoatomcut, rhocut, nonlocal potential cutoff
    % should be given by a table according to the element type
    rhoatomcut = 4.0;
    rhocut = 6.0;
    nlcut = 2.0;

    tempEntry.cutoffs.rad = rhocut;
    tempEntry.cutoffs.vLocal = rhocut;
    tempEntry.cutoffs.drv_vLocal = rhocut;
    tempEntry.cutoffs.pseudoCharge = rhocut;
    tempEntry.cutoffs.drv_pseudoCharge = rhocut;
    
    tempEntry.cutoffs.rhoAtom = rhoatomcut;
    tempEntry.cutoffs.drv_rhoAtom = rhoatomcut;

    for j = 1 : upf_nproj     
        tempEntry.cutoffs.nonlocal(j) = nlcut;
    end

    
    % ------------------- <PP_MESH> ------------------------------
    % read mesh
    read_start_field("PP_MESH", upfid);
    
    read_start_field("PP_R", upfid);
    buff_str = read_content_field("PP_R", upfid);
    upf_r = sscanf(buff_str, '%f', upf_mesh_size);
    
    read_start_field("PP_RAB", upfid);
    buff_str = read_content_field("PP_RAB", upfid);
    upf_rab = sscanf(buff_str, '%f', upf_mesh_size);
    
    read_end_field("PP_MESH", upfid);
    
    % add the mesh into samples
    tempEntry.samples.rad = upf_r;
    
    
    % ---------------------- <PP_NLCC> --------------------------------
    % <PP_NLCC> may not exist
    if upf_nlcc_flag == "T"
        read_start_field("PP_NLCC", upfid);
        buff_str = read_content_field("PP_NLCC", upfid);
        upf_nlcc = sscanf(buff_str, '%f', upf_mesh_size);
    end

    
    % --------------------- <PP_LOCAL> --------------------------------
    read_start_field("PP_LOCAL", upfid);
    buff_str = read_content_field("PP_LOCAL", upfid);
    upf_vloc = sscanf(buff_str, '%f', upf_mesh_size);
    
    % add vlocal into the samples
    % vlocal derivative is 0
    upf_vloc = 0.5 * upf_vloc;
    
    % interpolation
    [r, vr] = spline_rad(upf_r, upf_vloc, 1);
    [splb, splc, spld] = spline(r, vr);
    upf_vloc = seval(upf_r, r, vr, splb, splc, spld);
    
    tempEntry.samples.vLocal = upf_vloc;    
    tempEntry.samples.pseudoCharge = upf_vloc;
        
    
    % -------------------- <PP_NONLOCAL> ------------------------
    read_start_field("PP_NONLOCAL", upfid);
    upf_proj_l = zeros(upf_nproj, 1);
    
    for j = 1 : upf_nproj
        item_name = "PP_BETA." + num2str(j); 
        tag = read_start_field(item_name, upfid);
        
        % get attribute
        index = get_attr(tag, "index");
        InfoPrint(0, ' index = %s \n', index);
        index = str2double(index);
        angular_momentum = get_attr(tag, "angular_momentum");
        InfoPrint(0, ' angular_momentum = %s \n' ,angular_momentum);
        angular_momentum = str2double(angular_momentum);
        if angular_momentum > upf_lmax
            error('angular momentum mistake in UPF file');
        end
        upf_proj_l(index) = angular_momentum;
        
        % read non_local part data
        buff_str = read_content_field(item_name, upfid);
        upf_vnl = sscanf(buff_str, '%f', upf_mesh_size);
        
        if mod(angular_momentum, 2) == 0
            [r, vr] = spline_rad(upf_r, upf_vnl, 0);
        else
            [r, vr] = spline_rad(upf_r, upf_vnl, 1);
        end
        vr = vr ./ r;
        [splb, splc, spld] = spline(r, vr);
        upf_vnl = seval(upf_r, r, vr, splb, splc, spld);        

        % write nonlocal
        % nonlocal derivative is 0
        tempEntry.samples.nonlocal(:, 2*j-1) = upf_vnl;
    end
    
    % compute number of projectors for each l
    % nproj_l(l) is the number of the projectors having angular
    % momentum l
    nproj_l = zeros(upf_lmax+1, 1);
    for k = 0 : upf_lmax
        nproj_l(k+1) = sum(upf_proj_l == k);
    end

    % read <PP_DIJ>
    tag = read_start_field("PP_DIJ", upfid);
    size = get_attr(tag, "size");
    InfoPrint(0, ' PP_DIJ size = %s \n\n', size);
    size = str2double(size);
    
    if size ~= upf_nproj*upf_nproj
        error('Number of non-zeros Dij differs from number of projectors.');
    end
    
    upf_ndij = size;
    buff_str = read_content_field("PP_DIJ", upfid);
    upf_d = sscanf(buff_str, '%f', upf_ndij);
    
    read_end_field("PP_NONLOCAL", upfid);
    
    % add Dij to the weights
    upf_d = reshape(upf_d, upf_nproj, []);
    % check if Dij has non-diagonal element
    if ~isdiag(upf_d)
        error('Non-local Dij has off-diagonal elements');
    end
    for j = 1 : upf_nproj     
    % FIXME: nonlocal cutoffs should be given by a table according to
    % the element type
        tempEntry.NLweights(j) = 0.5 * upf_d(j, j);        
        tempEntry.NLtypes(j) = upf_proj_l(j);
        
        tempEntry.cutoffs.nonlocal(j) = nlcut;
    end
    
    % --------------------- <PP_RHOATOM> -----------------------------
    read_start_field("PP_RHOATOM", upfid);
    buff_str = read_content_field("PP_RHOATOM", upfid);
    upf_rho_atom = sscanf(buff_str, '%f', upf_mesh_size);
    
    % add the spline part
    [r, vr] = spline_rad(upf_r, upf_rho_atom, 1);
    vr = vr ./ (4 * pi * r .* r);
    [splb, splc, spld] = spline(r, vr);
    upf_rho_atom = seval(upf_r, r, vr, splb, splc, spld);
    
    % add rho_atom into the samples
    % rho_atom derivative is 0.0
    tempEntry.samples.rhoAtom = upf_rho_atom;
        
    fclose(upfid);
end

end


%*********************************************************************%
%                        Sub Functions                                %
%*********************************************************************%


% read the start block of field "name", i.e., <name ...>
function buff = read_start_field(name, upfid)
    search_str = '<' + name;
    token = fgetl(upfid);
    while ~contains(token, search_str) && ~feof(upfid)
        token = fgetl(upfid);
    end
    if feof(upfid)
        error('EOF reached before start of field %s \n', name);
    end
    buff = token;
    if buff(end) == '>'
        return;
    end
        
    % read until '>' is found
    while buff(end) ~= '>' && ~feof(upfid)
        token = fgetl(upfid);
        buff = strcat(buff, token);
    end
    if feof(upfid)
        error('EOF reached before > of field %s \n', name);
    end
end

% read string only concludes content of the field
% upfid starts next to <name>, ends at </name>
function buff_str = read_content_field(name, upfid)
    buff_str = "";
    token = fgetl(upfid);
    end_str = "</" + name + ">"; 
    while ~contains(token, end_str)
        buff_str = buff_str + token;
        token = fgetl(upfid);
    end
end

% read the end block of field "name" if there exists, i.e., </name>
function read_end_field(name, upfid)
    search_str = "</" + name + ">";
    token = fgetl(upfid);
    while ~contains(token, search_str) && ~feof(upfid)
        token = fetl(upfid);
    end
    if feof(upfid)
        error('EOF reached before end of field %s \n', name);
    end
end

% search value of attribute "attr" in buff
function value = get_attr(buff, attr)
    search_str = " " + attr + "=";
    
    buff = char(buff);
    % find attribute name in buff
    idx_attr = strfind(buff, search_str);    
    if isempty(idx_attr)
        error('get_attr: attribute not found: %s \n', attr);
    else
        idx_value = strfind(buff(idx_attr : end), '"');
        if length(idx_value) < 2
            error('get_attr: attribute not found: %s \n', attr);
        end
        idx_begin = idx_attr - 1 + idx_value(1);
        idx_end = idx_attr  - 1 + idx_value(2);

        value = buff(idx_begin+1 : idx_end-1);
    end
end
    
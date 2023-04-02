function [tempEntry, atomicNum] = ReadUPF(fileName)
% READUPF reads pseudopotential info from UPF file
%
%    [tempEntry, atomicNum] = ReadUPF(fileName) reads data from file and
%    returns the atomic number, pseudopotential data stored in the 
%    structure PTEntry.
%
%    See also PeriodTable, PTEntry.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

% NOTE:
% Currently only version 2 UPF file with norm-conserving condition can be
% supported.

etable = ElementTable();
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

    % -------------------- <PP_HEADER> ------------------------------
    tag = read_start_field("PP_HEADER", upfid);
    
    % get attribure "element"
    upf_symbol = get_attr(tag, "element");
    upf_symbol = erase(upf_symbol, ' ');
    
    % get atomic number and mass
    atomicNum = etable.z(upf_symbol);
    mass = etable.mass(upf_symbol);  
    
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
    upf_zval = str2double(upf_zval);
        
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


    %
    % Initialize tempEntry
    %
    tempEntry.params.Znuc = atomicNum;
    tempEntry.params.Mass = mass;
    tempEntry.params.Zval = upf_zval;
    
    tempEntry.samples.nonlocal = zeros(upf_mesh_size, 2*upf_nproj);
    tempEntry.cutoffs.nonlocal = zeros(upf_nproj, 1);
    tempEntry.nlweights        = zeros(upf_nproj, 1);
    tempEntry.nltypes          = zeros(upf_nproj, 1);
    
    tempEntry.samples.drv_vLocalSR     = zeros(upf_mesh_size, 1);
    tempEntry.samples.drv_rhoAtom      = zeros(upf_mesh_size, 1);    
    
    rhocut = get_attr(tag, "rho_cutoff");
    rhocut = str2double(rhocut);
    tempEntry.cutoffs.rad = rhocut;
    tempEntry.cutoffs.vLocalSR = rhocut;
    tempEntry.cutoffs.drv_vLocalSR = rhocut;
    tempEntry.cutoffs.rhoAtom = rhocut;
    tempEntry.cutoffs.drv_rhoAtom = rhocut;
    for j = 1 : upf_nproj     
        tempEntry.cutoffs.nonlocal(j) = rhocut;
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

    % if rho_cutoff is not provided in UPF file
    if rhocut == 0
        rhocut = upf_r(end);
        tempEntry.cutoffs.rad = rhocut;
        tempEntry.cutoffs.rhoAtom = rhocut;
        tempEntry.cutoffs.drv_rhoAtom = rhocut;
        for j = 1 : upf_nproj     
            tempEntry.cutoffs.nonlocal(j) = rhocut;
        end    
    end
    
    
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
    upf_vloc = 0.5 * upf_vloc;  % unit conversion: from Ry to Bohr
    % interpolation
    [r, vr] = spline_rad(upf_r, upf_vloc, 1);
    [splb, splc, spld] = spline(r, vr);
    upf_vloc = seval(upf_r, r, vr, splb, splc, spld);
    tempEntry.samples.vLocalSR = upf_vloc;    
    % NOTE: here is still initial data, a post-processing in 
    % PeriodTable/Setup() will generate real short range part

    %
    % determine RGaussian through an optimization problem
    %
    fitfunc = fittype('-Zval*erf(r/Rg)/r', 'problem', 'Zval', 'independent', 'r');
    RGaussian_threshold = 1e-4;
    % find the start point to perform curve fitting
    idx_position = floor( length(upf_r) * 5/6 );
    idx_threshold = ...
        find( abs(upf_vloc - (-upf_zval./upf_r)) <= RGaussian_threshold );
    if isempty(idx_threshold)
        idx = idx_position; 
    else
        % skip the potential oscillation
        idx_increment = diff(idx_threshold);
        temp = find(idx_increment > 1);
        if ~isempty(temp)
            idx_threshold = idx_threshold((temp(end)+1) : end);
        end
        idx = idx_threshold;
    end
    fitParam = fit(upf_r(idx : end), upf_vloc(idx : end), fitfunc, ...
                    'problem', upf_zval, 'startpoint', 1.0, 'lower', 0.0);
    RGaussian = fitParam.Rg;
    tempEntry.params.RGaussian = RGaussian;
    
    %
    % determine vsrcut according to RGaussian
    %
    vlocSR_threshold = 1e-6;  % threshold for short range part of VLocal
    vGaussian = -upf_zval * erf(upf_r ./ RGaussian) ./ upf_r;
    idx_threshold = find( abs(upf_vloc - vGaussian) <= vlocSR_threshold );
    if isempty(idx_threshold)
        vsrcut = upf_r(end);
    else
        % skip the potential osillation
        idx_increment = diff(idx_threshold);
        temp = find(idx_increment > 1);
        if ~isempty(temp)
            idx_threshold = idx_threshold((temp(end)+1) : end);
        end
        vsrcut = upf_r(idx_threshold(1));
    end
    tempEntry.cutoffs.vLocalSR = vsrcut;
    tempEntry.cutoffs.drv_vLocalSR = vsrcut;

    InfoPrint(0, ' RGaussian = ', tempEntry.params.RGaussian);
    InfoPrint(0, ' vsrcut = ' , tempEntry.cutoffs.vLocalSR);

    
    % -------------------- <PP_NONLOCAL> ------------------------
    read_start_field("PP_NONLOCAL", upfid);
    upf_proj_l = zeros(upf_nproj, 1);
    
    % read <PP_BETA>
    for j = 1 : upf_nproj
        item_name = "PP_BETA." + num2str(j); 
        tag = read_start_field(item_name, upfid);
        
        % get attribute
        index = get_attr(tag, "index");
        InfoPrint(0, ' index = %s \n', index);
        index = str2double(index);
        angular_momentum = get_attr(tag, "angular_momentum");
        InfoPrint(0, ' angular_momentum = %s \n', angular_momentum);
        angular_momentum = str2double(angular_momentum);
        if angular_momentum > upf_lmax
            error('angular momentum mistake in UPF file');
        end
        upf_proj_l(index) = angular_momentum;
        
        cutoff_radius = get_attr(tag, "cutoff_radius");
        InfoPrint(0, ' cutoff_radius = %s \n', cutoff_radius);
        tempEntry.cutoffs.nonlocal(j) = str2double(cutoff_radius);
        
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

    tempEntry.nltypes = upf_proj_l;
    
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
        
    % add Dij to the weights
    upf_d = reshape(upf_d, upf_nproj, []);
    upf_d = 0.5 * upf_d;  % unit conversion from Ry to Bohr
    tempEntry.nlweights = diag(upf_d);
    
    % Post processing nonlocal pseudopotential by diagonalizing
    % nonlocal D and combine the vectors into upf_vnl.
    % This is done by diagonalizing each angular momentum block.
    for l = 0 : upf_lmax
        % extract relevant index
        idx = find(upf_proj_l == l);
        D_block = upf_d(idx, idx);
        % If D_block is a diagonal matrix, no need for diagonalization,
        % just copy the vectors to right place.
        % If D_block has non-diagonal element, diagonalize the sublock
        % and combine eigenvectors to projectors.
        if ~isdiag(D_block)
            [tempV, tempD] = eig(D_block);
            vnl_block = tempEntry.samples.nonlocal(:, 2*idx-1);
            vnl_block = vnl_block * tempV;
            tempEntry.samples.nonlocal(:, 2*idx-1) = vnl_block;
            tempEntry.nlweights(idx) = diag(tempD);
        end    
    end

    read_end_field("PP_NONLOCAL", upfid);
    
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
    token = fgets(upfid);
    end_str = "</" + name + ">"; 
    while ~contains(token, end_str)
        buff_str = buff_str + token;
        token = fgets(upfid);
    end
end

% read the end block of field "name" if there exists, i.e., </name>
% NOTE: this function is used to read end line of field that contains 
% subfields, e.g., PP_NONLOCAL. Fields without subfield should NOT 
% use this function to read end line of the field.
function read_end_field(name, upfid)
    search_str = "</" + name + ">";
    token = fgetl(upfid);
    while ~contains(token, search_str) && ~feof(upfid)
        token = fgetl(upfid);
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
        % add another try to support different types of UPF file
        search_str = attr + "=";
        idx_attr = strfind(buff, search_str);
    end    
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

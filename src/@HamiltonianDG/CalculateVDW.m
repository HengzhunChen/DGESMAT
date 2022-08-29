function [EVdw, forceVdw] = CalculateVDW(HamDG, VDWType)
% HAMILTONIANDG/CALCULATEVDW calculates van der Waals energy
%    and force.
%
%    [EVdw, forceVdw] = CalculateVDW(HamDG, VDWType) computes the van der
%    Waals energy and force according to the type of VDW VDWType.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

% NOTE: only some types of exchange-correlation are supported here.

% TODO: add support for more types of exchange-correlation

atomList = HamDG.atomList;
dm = HamDG.domain;
EVdw = 0;
forceVdw = zeros(length(atomList), dimDef());

if VDWType == "DFT-D2"
    % vdw_nspecies = 55;
    vdw_d = 20.0;
    vdw_tol_default = 1e-10;
    vdw_s_pbe = 0.75;
    vdw_s_blyb = 1.2;
    vdw_s_b3lyp = 1.05;
    vdw_s_hse = 0.75;
    vdw_s_pbe0 = 0.60;
    % Thin Solid Films 535 (2013) 387-389
    % J. Chem. Theory Comput. 2011, 7, 88–96
    
    vdw_c6_dftd2 = ...
    [ 0.14, 0.08, 1.61, 1.61, 3.13, 1.75, 1.23, 0.70, 0.75, 0.63, ...
      5.71, 5.71,10.79, 9.23, 7.84, 5.57, 5.07, 4.61,10.80,10.80, ...
      10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80, ...
      16.99,17.10,16.37,12.64,12.47,12.01,24.67,24.67,24.67,24.67, ...
      24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,37.32,38.71, ...
      38.44,31.74,31.50,29.99, 0.00 ];
  
    vdw_r0_dftd2 = ...
    [ 1.001,1.012,0.825,1.408,1.485,1.452,1.397,1.342,1.287,1.243, ...
      1.144,1.364,1.639,1.716,1.705,1.683,1.639,1.595,1.485,1.474, ...
      1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562, ...
      1.650,1.727,1.760,1.771,1.749,1.727,1.628,1.606,1.639,1.639, ...
      1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.672,1.804, ...
      1.881,1.892,1.892,1.881,1.000 ];

    vdw_c6_dftd2 = vdw_c6_dftd2 ./ 2625499.62 .* ((10/0.52917706) ^6);
    vdw_r0_dftd2 = vdw_r0_dftd2 ./ 0.52917706;
    
    %   vdw_c6(i, j) = sqrt( vdw_c6_dftd2(i) * vdw_c6_dftd2(j) );
    %   vdw_r0(i, j) = vdw_r0_dftd2(i) + vdw_r0_dftd2(j);
    vdw_c6 = sqrt(vdw_c6_dftd2' * vdw_c6_dftd2);
    vdw_r0 = vdw_r0_dftd2' + vdw_r0_dftd2;
    
    if HamDG.XCType == "XC_GGA_XC_PBE"
      vdw_s = vdw_s_pbe;
    else
      error("Van der Waals DFT-D2 correction in only compatible with GGA-PBE!");
    end
    
    for ii = -1 : 1
        for jj = -1 : 1
            for kk = -1 : 1
                
                for i = 1 : length(atomList)
                    iType = atomList(i).type;
                    for j = 1 : i
                        jType = atomList(j).type;
                        
                        rx = atomList(i).pos(1) - atomList(j).pos(1) + ii * dm.length(1);
                        ry = atomList(i).pos(2) - atomList(j).pos(2) + jj * dm.length(2);
                        rz = atomList(i).pos(3) - atomList(j).pos(3) + kk * dm.length(3);
                        rr = sqrt(rx * rx + ry * ry + rz * rz);
                        
                        if rr > 0.0001 && rr < 75.0
                            sfactor = vdw_s;
                            if i == j
                                sfactor = sfactor * 0.5;
                            end
                            
                            c6 = vdw_c6(iType, jType);
                            r0 = vdw_r0(iType, jType);
                            
                            ex = exp( -vdw_d * (rr / r0 - 1) );
                            fr = 1 / (1 + ex);
                            c6r6 = c6 / rr^6;
                            
                            % Contribution to energy
                            EVdw = EVdw - sfactor * fr * c6r6;
                            
                            % Contribution to force
                            if i ~= j
                            
                                gr = (vdw_d / r0) * (fr * fr) * ex;
                                grad = sfactor * (gr - 6.0 * fr / rr) * c6r6 / rr;
                                
                                fx = grad * rx;
                                fy = grad * ry;
                                fz = grad * rz;
                                
                                forceVdw(i, 1) = forceVdw(i, 1) + fx;
                                forceVdw(i, 2) = forceVdw(i, 2) + fy;
                                forceVdw(i, 3) = forceVdw(i, 3) + fz;
                                forceVdw(j, 1) = forceVdw(j, 1) - fx;
                                forceVdw(j, 2) = forceVdw(j, 2) - fy;
                                forceVdw(j, 3) = forceVdw(j, 3) - fz;
                            
                            end  % end for i ~= j
                        end  % end if
                        
                    end  % end for j
                end  % end for i
                
            end  % end for ii
        end  % end for jj
    end   % end for kk
    
end  % end if DFT-D2

end
function [Eself, EIonSR, forceIonSR] = CalculateIonSelfEnergyAndForce(ptable, atomList, domain)
% PERIODTABLE/CALCULATEIONSELFENERGYANDFORCE calculates self energy and
%    short range ionic energy and force.
%
%    [Eself, EIonSR, forceIonSR] = CalculateIonSelfEnergyAndForce(ptable, atomList, domain)
%    computes the self energy and short range ionic energy and force
%    according to atoms in atomList over domain.  
%
%    See also Domain, Atom, HamiltonianKS/CalculatePseudoPotential, 
%    HamiltonianDG/CalculatePseudoPotential.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numAtom = length(atomList);

% ----------------------- Self energy ----------------------------------
Eself = 0;
for i = 1 : numAtom
    type = atomList(i).type; 
    Eself = Eself + ptable.SelfIonInteraction(type);
end

% ------------------ Short range repulsion part ------------------------
EIonSR = 0;
forceIonSR = zeros(numAtom, dimDef());

dm = domain;

for a = 1 : numAtom
    type_a = atomList(a).type;
    Zval_a = ptable.Zval(type_a);
    RGaussian_a = ptable.RGaussian(type_a);
    
    for b = a : numAtom
        same_atom = a == b;
        % need to consider the interaction between the same atom and
        % its periodic image. Be sure not to double count
        
        type_b = atomList(b).type;
        Zval_b = ptable.Zval(type_b);
        RGaussian_b = ptable.RGaussian(type_b);
        
        radius_ab = sqrt( RGaussian_a*RGaussian_a + RGaussian_b*RGaussian_b );
        
        % convergence criterion for lattice sums:
        % facNbr * radius < ncell * d
        
        facNbr = 8.0;
        ncell1 = floor( facNbr * radius_ab / dm.length(1) );
        ncell2 = floor( facNbr * radius_ab / dm.length(2) );
        ncell3 = floor( facNbr * radius_ab / dm.length(3) );
        
        pos_ab = atomList(a).pos - atomList(b).pos;
        pos_ab = pos_ab - round(pos_ab ./ dm.length) .* dm.length;
        
        % loop over neighboring cells
        for ic1 = -ncell1 : ncell1
            for ic2 = -ncell2 : ncell2
                for ic3 = -ncell3 : ncell3
                    if ~same_atom || ic1 ~= 0 || ic2 ~= 0 || ic3 ~= 0
                        if same_atom
                            fac = 0.5;
                        else
                            fac = 1.0;
                        end
                        
                        pos_ab_image(1) = pos_ab(1) + ic1 * dm.length(1);
                        pos_ab_image(2) = pos_ab(2) + ic2 * dm.length(2);
                        pos_ab_image(3) = pos_ab(3) + ic3 * dm.length(3);
                        
                        r_ab = norm(pos_ab_image, 2);
                        esr_term = Zval_a * Zval_b * erfc(r_ab / radius_ab) / r_ab;
                        desr_erfc = 2 * Zval_a * Zval_b * ...
                            exp( -(r_ab / radius_ab)^2 ) / (radius_ab * sqrt(pi));
                        % desrdr = (1/r) d Esr / dr
                        desrdr = - fac * (esr_term + desr_erfc) / (r_ab*r_ab);
                        
                        EIonSR = EIonSR + fac * esr_term;
                        
                        forceIonSR(a, 1) = forceIonSR(a, 1) - desrdr * pos_ab_image(1);
                        forceIonSR(a, 2) = forceIonSR(a, 2) - desrdr * pos_ab_image(2);
                        forceIonSR(a, 3) = forceIonSR(a, 3) - desrdr * pos_ab_image(3);
                        forceIonSR(b, 1) = forceIonSR(b, 1) + desrdr * pos_ab_image(1);
                        forceIonSR(b, 2) = forceIonSR(b, 2) + desrdr * pos_ab_image(2);
                        forceIonSR(b, 3) = forceIonSR(b, 3) + desrdr * pos_ab_image(3);
                    end
                    
                end  % end for ic1
            end  % end for ic2
        end  % end for ic3
        
    end  % end for b
end  % end for a
    
end
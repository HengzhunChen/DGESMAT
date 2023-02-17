function [occupationRate, fermi] = CalculateOccupationRate(scf, eigVal)
% SCF/CALCULATEOCCUPATIONRATE For a given finite temperature, update the 
%    occupation number and Fermi energy according to eigenvalue eigVal.
%
%    NOTE: if number of occupied state equals to total number of states,
%    occupationRate is an all one column vector, else it will have the 
%    same size as eigVal.
%
%    See also SCF.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% Magic number here
tol = 1e-10;
maxiter = 100;

hamKS = scf.eigSol.hamKS;
npsi = hamKS.NumStateTotal();
nOccStates = hamKS.numOccupiedState;
Tbeta = scf.Tbeta;

neig = length(eigVal);
if neig ~= npsi
    message = "The number of eigensates do not match. " + ...
              "num of eigVal = " + num2str(neig) + ...
              " numStateTotal  = " + num2str(npsi);
    error(message);
end


if npsi > nOccStates
    % use bisection to find efermi such that 
    %
    % sum_i fermi_dirac( ev(i) ) = nooc
    %
    
    ilb = nOccStates - 1;
    iub = nOccStates + 1;
    if ilb <= 0
        message = ...
            "The chemical potential is smaller than the lowest eigvalue. " + ...
            "Please set Extra_States = 0 to avoid this bug.";
        error(message);
    end
    lb = eigVal(ilb);
    ub = eigVal(iub);
    
    % Calculate Fermi-Dirac function and make sure that 
    % flb < nooc and fub > nocc
    
    flb = sum( 1 ./ (1 + exp( Tbeta*(eigVal - lb) )) );
    fub = sum( 1 ./ (1 + exp( Tbeta*(eigVal - ub) )) );
    
    while (nOccStates-flb) * (fub-nOccStates) < 0
        if flb >  nOccStates
            if ilb > 1
                ilb = ilb - 1;
                lb = eigVal(ilb);
                flb = sum( 1 ./ (1 + exp(Tbeta*(eigVal - lb) )) );
            else
                error('Cannot find a lower bound for efermi');
            end
        end
        if fub < nOccStates
            if iub < npsi-1
                iub = iub + 1;
                ub = eigVal(iub);
                fub = sum( 1 ./ (1 + exp(Tbeta*(eigVal - ub) )) );
            else
                error('Cannot find a lower bound for efermi, try to increase the number of wavefunctions');
            end
        end
    end
    
    fermi = (lb + ub) * 0.5;
    occupationRate = 1 ./ (1 + exp(Tbeta*(eigVal - fermi)));
    occsum = sum(occupationRate);
    
    % start bisection iteration
    iter = 1;
    while abs(occsum - nOccStates) > tol && iter < maxiter
        if occsum < nOccStates
            lb = fermi;
        else
            ub = fermi;
        end
        
        fermi = (lb + ub) * 0.5;
        occupationRate = 1 ./ (1 + exp(Tbeta*(eigVal - fermi)));
        occsum = sum(occupationRate);
        
        iter = iter + 1;
    end
    
else
    if npsi == nOccStates
        occupationRate = ones(npsi, 1);
        fermi = eigVal(npsi);
    else
        error('The number of eigenvalues in eigVal should be larger than nocc');
    end
end

end
function [occupationRate, fermi] = CalculateOccupationRate(scfDG, eigVal)
% SCFDG/CALCULATEOCCUPATIONRATE For a given finite temperature, update the 
%    occupation number and Fermi energy according to eigenvalue eigVal.
%
%    NOTE: if number of occupied state equals to total number of states,
%    occupationRate is an all one column vector, else it will have the 
%    same size as eigVal.
%
%    See also SCFDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% TODO: the case when npsi > nOccStates have NOT be tested

% FIXME Magic number here
tol = 1e-10;
maxiter = 100;

Tbeta = scfDG.Tbeta;
Tsigma = scfDG.Tsigma;

SmearingScheme = scfDG.smearing.SmearingScheme;
MPsmearingOrder = scfDG.smearing.MPsmearingOrder;

if SmearingScheme == "FD"
    npsi = scfDG.hamDG.NumStateTotal();
    nOccStates = scfDG.hamDG.numOccupiedState;
    
    neig = length(eigVal);
    if neig ~= npsi
        message = "The number of eigenstates do not match. " + ...
                  "eigVal = " + num2str(neig) + ...
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
            error('The number of eigenvalues in ev should be larger than nocc');
        end    
    end
    % end of SmearingScheme == "FD"
else
    % MP and GB type smearing
    npsi = scfDG.hamDG.NumStateTotal();
    nOccStates = scfDG.hamDG.numOccupiedState;
    
    neig = length(eigVal);
    if neig ~= npsi
        message = "The number of eigensates do not match. " + ...
                  "eigVal = " + num2str(neig) + ...
                  " numStateTotal  = " + num2str(npsi);
        error(message);
    end

    if npsi > nOccStates
        % set up the bounds
        lb = eigVal(1);
        ub = eigVal(end);
        
        % set up the function bounds
        flb = MPoccupResidual(Tsigma, MPsmearingOrder, eigVal, lb, nOccStates);
        fub = MPoccupResidual(Tsigma, MPsmearingOrder, eigVal, ub, nOccStates);
        
        if flb * fub > 0.0
            error('Bisection method for finding Fermi level cannot proceed !!');
        end
        
        fermi = (lb + ub) * 0.5;
        
        % start bisection iteration
        iter = 1;
        fx = MPoccupResidual(Tsigma, MPsmearingOrder, eigVal, fermi, nOccStates);
        
        % iterate using the bisection method
        while abs(fx) > tol && iter < maxiter
            flb = MPoccupResidual(Tsigma, MPsmearingOrder, eigVal, lb, nOccStates);
            
            if flb * fx < 0.0
                ub = fermi;
            else
                lb = fermi;
            end
            
            fermi = (lb + ub) * 0.5;
            fx = MPoccupResidual(Tsigma, MPsmearingOrder, eigVal, fermi, nOccStates);
            
            iter = iter + 1;
        end
        
        if iter >= maxiter
            error('Bisection method for finding Fermi level does not appear to converge !!');
        else
            % Bisection method seem to have converged
            % Fill up the occupations
            occupationRate = PopulateMPoccup(Tsigma, MPsmearingOrder, eigVal, fermi);
        end
        % end of if (npsi > nOccStates)
    else
        if npsi == nOccStates
            occupationRate = ones(npsi, 1);
            fermi = eigVal(npsi);
        else
            % npsi < nOccStates
            error('The number of top eigenvalues should be larger than number of occupied states !!');
        end
    end
end    

end
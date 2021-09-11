%% Heat transfer correlation by Gelperin / Einstein
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.; Schwarzmayr, P.  
%Experimental Investigation of the Heat Transfer between Finned Tubes and 
%a Bubbling Fluidized Bed with Horizontal Sand Mass Flow. Energies 2021, 
%14, x. https://doi.org/10.3390/xxxxx
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5500329
%
%
%
%This function describes the heat transfer correlation in a fluidized bed
%according to:
%
%Gelperin, N.I.; Einstein V.G. Heat Transfer in Fluidized Beds. In 
%Fluidization; Davidson, J.F., Harrison, D., Eds.; Academic Press: London, 
%United Kingdom, 1971; pp. 471-540. https://doi.org/10.1002/aic.690330123
%
%
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir


function alpha=gelperinEinstein(d_p,rho_p,d_tube,pitch_horz,pitch_vert,T)
    persistent g 
    if isempty(g)
        g=9.81;
    end
    
    lambda_g=DryAir.lambda(T);
    rho_g=DryAir.rho(T);
    eta_g=DryAir.eta(T);
    
    Ar=rho_g.*d_p.^3.*(rho_p-rho_g).*g./eta_g.^2;
    
    Nu_p=0.74.*Ar.^0.22.*(1-d_tube./pitch_vert.*(1+d_tube./(pitch_horz+d_tube))).^0.25;
    
    alpha=Nu_p.*lambda_g./d_p;
end





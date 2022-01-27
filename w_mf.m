%% Minimum fluidization velocity function
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
%This function describes the minimum fluidization velocity according to
%Richardson:
%
%Kunii, D.; Levenspiel, O. Fluidization and Mapping of Regimes. In 
%Fluidization Engineering, 2nd ed.; Brenner, H., Ed.; 
%Butter-worth-Heinemann: Boston, Massachusetts, United States, 1991; pp. 
%61-94. https://doi.org/10.1016/B978-0-08-050664-7.50009-3
%
%
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir


function w_mf=w_mf(d_p,rho_p,p,T)
    persistent C1 C2 g
    if isempty(C1)
        %C1, C2 according to Richardson
        C1=25.7;
        C2=0.0365;
        
        g=9.81;
    end
    
    
    eta_g=DryAir.eta(T);
    rho_g=DryAir.rho(p,T);
    
    Ar=rho_g.*d_p.^3.*(rho_p-rho_g).*g./eta_g.^2;
    Re=sqrt(C1.^2+C2.*Ar)-C1;
    
    w_mf=Re.*eta_g./(d_p.*rho_g);
end





%% Bed porosity function
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
%https://doi.org/10.5281/zenodo.5802407
%
%
%
%This function describes the porosity in a fluidized bed according to
%Goroshko:
%
%Zabrodsky, S.S. Hydrodynamics and Heat Transfer in Fluidized Beds, 
%7th ed.; MIT Press: Cambridge, Massachusetts, United States, 1966; p. 71.
%
%
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir
%   - w_mf.m


function eps=epsilon(d_p,rho_p,FG,p_A,T_A)
    persistent g 
    if isempty(g)
        g=9.81;
    end
    
    eta_A=DryAir.eta(T_A);
    rho_A=DryAir.rho(p_A,T_A);
    
    wmf=w_mf(d_p,rho_p,p_A,T_A);
    w=wmf.*FG;
    
    Ar=rho_A.*d_p.^3.*(rho_p-rho_A).*g./eta_A.^2;
    Re=d_p.*rho_A.*w./eta_A;
    
    eps=((18.*Re+0.36.*Re.^2)./Ar).^0.21;
end





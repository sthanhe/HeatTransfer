%% Heat transfer correlation by Molerus
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
%Molerus, O.; Burschka, A.; Dietz, S. Particle migration at solid surfaces 
%and heat transfer in bubbling fluidized beds—II. Prediction of heat 
%transfer in bubbling fluidized beds. Chem. Eng. Sci. 1995, 50, 879-885. 
%https://doi.org/10.1016/0009-2509(94)00446-X
%
%
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
%   - @SiO2


function alpha=molerus(rho_p,p_A,T_A)
    persistent g
    if isempty(g)
        g=9.81;
    end
    
    lambda_A=DryAir.lambda(T_A);
    eta_A=DryAir.eta(T_A);
    rho_A=DryAir.rho(p_A,T_A);
    
    l=(eta_A./(sqrt(g).*(rho_p-rho_A))).^(2/3);
    alpha=(l./(0.09.*lambda_A)+l./(0.09.*SiO2.c_p(T_A).*eta_A)).^-1;
end





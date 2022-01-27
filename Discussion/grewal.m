%% Heat transfer correlation by Grewal
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
%Grewal, N.S. A generalized correlation for heat transfer between a
%gas—solid fluidized bed of small particles and an im-mersed staggered 
%array of horizontal tubes. Powder Technol. 1981, 30, 145-154. 
%https://doi.org/10.1016/0032-5910(81)80007-1
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
%   - epsilon.m
%   - w_mf.m


function alpha=grewal(d_p,rho_p,FG,d_tube,pitch,p_A,T_A,eps)
    persistent g 
    if isempty(g)
        g=9.81;
    end
    
    if nargin<8
        eps=epsilon(d_p,rho_p,FG,p_A,T_A);
    end
    
    lambda_A=DryAir.lambda(T_A);
    cp_A=DryAir.c_p(T_A);
    eta_A=DryAir.eta(T_A);
    
    wmf=w_mf(d_p,rho_p,p_A,T_A);
    w=wmf.*FG;
    
    Nu_T=47.*(1-eps);
    Nu_T=Nu_T.*(w.*d_tube.*eta_A./(d_p.^3.*rho_p.*g)).^0.325;
    Nu_T=Nu_T.*(eta_A.*cp_A./lambda_A).^0.3;
    Nu_T=Nu_T.*(rho_p.*SiO2.c_p(T_A).*d_tube.^1.5*sqrt(g)./lambda_A).^0.23;
    Nu_T=Nu_T.*(1-0.21.*(pitch./d_tube).^-1.75);
    
    alpha=Nu_T.*lambda_A./d_tube;
end





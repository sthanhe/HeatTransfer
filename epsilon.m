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
%https://doi.org/10.5281/zenodo.5500329
%
%
%
%This function describes the porosity in a fluidized bed required by the
%heat transfer correlation functions by Andeen / Glicksman and Grewal:
%
%Andeen, B.R.; Glicksman, L.R. Paper 76-HT-67. ASME-AIChE Heat Transfer 
%Conf.., St. Louis, Missouri, United States, Aug 9-11, 1976.
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
%Necessary files, classes and functions:
%   - @DryAir
%   - w_mf.m


function eps=epsilon(d_p,rho_p,FG,T)
    persistent g 
    if isempty(g)
        g=9.81;
    end
    
    rho_g=DryAir.rho(T);
    eta_g=DryAir.eta(T);
    
    wmf=w_mf(d_p,rho_p,T);
    w=wmf.*FG;
    
    Ar=rho_g.*d_p.^3.*(rho_p-rho_g).*g./eta_g.^2;
    Re_p=d_p.*rho_g.*w./eta_g;
    
    eps=((18.*Re_p+0.36.*Re_p.^2)./Ar).^0.21;
end





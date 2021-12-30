%% Heat transfer correlation by Andeen / Glicksman
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
%This function describes the heat transfer correlation in a fluidized bed
%according to:
%
%Andeen, B.R.; Glicksman, L.R. Paper 76-HT-67. ASME-AIChE Heat Transfer 
%Conf.., St. Louis, Missouri, United States, Aug 9-11, 1976.
%
%
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir
%   - epsilon.m
%   - w_mf.m


function [alpha,Nu_T]=andeenGlicksman(d_p,rho_p,FG,d_tube,p_A,T_A,eps)
    persistent g 
    if isempty(g)
        g=9.81;
    end
    
    if nargin<7
        eps=epsilon(d_p,rho_p,FG,p_A,T_A);
    end
    
    lambda_A=DryAir.lambda(T_A);
    cp_A=DryAir.c_p(T_A);
    eta_A=DryAir.eta(T_A);
    
    wmf=w_mf(d_p,rho_p,p_A,T_A);
    w=wmf.*FG;
    
    Nu_T=900.*(1-eps);
    Nu_T=Nu_T.*(w.*d_tube.*eta_A./(d_p.^3.*rho_p.*g)).^0.326;
    Nu_T=Nu_T.*(eta_A.*cp_A./lambda_A).^0.3;
    
    alpha=Nu_T.*lambda_A./d_tube;
end





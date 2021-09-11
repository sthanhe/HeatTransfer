%% Heat transfer correlation by Zabrodsky
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
%Zabrodsky, S.S.; Antonishin, N.V.; Parnas, A.L. On fluidized 
%bed-to-surface heat transfer. Can. J. Chem. Eng. 1976, 54, 52-58. 
%https://doi.org/10.1002/cjce.5450540107
%
%
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir


function alpha=zabrodsky(d_p,rho_p,T)
    alpha=35.8.*rho_p.^0.2.*DryAir.lambda(T).^0.6.*d_p.^-0.36;
end





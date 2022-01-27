%% Heat transfer correlation by Martin
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
%Martin, H. Heat Transfer in Fluidized Beds. In VDI Heat Atlas, 2nd ed.;
%Stephan, P., Kabelac, S.., et al., Eds.; Springer: Berlin Heidelberg,
%Germany, 2010; pp. 1301–1310. https://doi.org/10.1007/978-3-540-77877-6_98
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


function [alpha,alpha_P,alpha_G,alpha_R]=martin(d_p,rho_p,lambda_p,FG,p_A,T_A,eps_mf,eps_R,eps)
    persistent C_A C_K g sigma
    if isempty(C_A)
        C_A=2.8;    %Air
        C_K=2.6;
        
        
        g=9.81;                 %gravitational acceleration
        sigma=5.670374419e-8;   %Stefan Boltzmann constant
    end
    
    if nargin<9
        eps=epsilon(d_p,rho_p,FG,p_A,T_A);
        eps(eps<eps_mf)=eps_mf;
    end
    
    lambda_A=DryAir.lambda(T_A);
    eta_A=DryAir.eta(T_A);
    cp_A=DryAir.c_p(T_A);
    Pr_A=DryAir.Pr(T_A);
    rho_A=DryAir.rho(p_A,T_A);
    
    gamma=NaN(size(T_A));
    for i=1:numel(T_A)
        gamma(i)=fzero(@(gamma) log10(1/gamma-1)-0.6+(1000/T_A(i)+1)/C_A,[1e-3,1-1e-3]);
    end
    Lambda=sqrt(2*pi*DryAir.R.*T_A).*lambda_A./(p_A.*(2*cp_A-DryAir.R));
    l=2*Lambda.*(2./gamma-1);
    
    Nu_WPmax=4*((1+2*l./d_p).*log(1+d_p./(2*l))-1);
    Z=rho_p.*SiO2.c_p(T_A)./(6*lambda_A).*sqrt(g.*d_p.^3.*(eps-eps_mf)./(5*(1-eps_mf).*(1-eps)));
    Nu_WP=(Nu_WPmax.^-1+lambda_A./lambda_p./(4*(1+sqrt(3*C_K*lambda_A.*Z./(2*pi*lambda_p))))).^-1;
    
    N=Nu_WP./(C_K.*Z);
    Nu_p=(1-eps).*Z.*(1-exp(-N));
    alpha_P=Nu_p.*lambda_A./d_p;
    
    
    Ar=rho_A.*(rho_p-rho_A)*g.*d_p.^3./(eta_A.^2);
    Nu_g=0.009*Pr_A.^(1/3).*sqrt(Ar);
    alpha_G=Nu_g.*lambda_A./d_p;
    
    
    alpha_R=4*eps_R*sigma*T_A.^3;
    
    
    alpha=alpha_P+alpha_G+alpha_R;
end





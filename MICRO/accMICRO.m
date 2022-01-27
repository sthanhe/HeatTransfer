%% Accuracy of MICRO test rig results
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
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.5890230
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5500329
%
%
%
%This function calculates the measurement uncertainties of all MICRO test 
%rig results.
%This function is only called during the general analysis of the results in
%the script "Analyze_MICRO.m".
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir
%   - w_mf.m


function tab=accMICRO(tab,d_p)
    %% Constants
    l=0.2;
    w=0.2;
    
    l_heated=0.1;
    d=25e-3;
    A=d*pi*l_heated;
    
    rho_p=2650;
    rho_Anorm=1.293;
    
    p_amb=101325;
    p_A=p_amb+4335;
    
    
    %% Calculate ranges
    [T_BED1min,T_BED1max]=TRange(tab.T_BED1);
    [T_PROBE1min,T_PROBE1max]=TRange(tab.T_PROBE1);
    [T_PROBE2min,T_PROBE2max]=TRange(tab.T_PROBE2);
    [T_LeinMin,T_LeinMax]=TRange(tab.T_Lein);
    [T_LausMin,T_LausMax]=TRange(tab.T_Laus);

    [VdotAir15min,VdotAir15max]=VdotRange(tab.VdotAir15,tab.VdotAir);
    [VdotAir20min,VdotAir20max]=VdotRange(tab.VdotAir20,tab.VdotAir);

    [Umin,Umax]=URange(tab.U);
    [Imin,Imax]=IRange(tab.I);


    %% Degree of Fluidization
    rho_Amin=DryAir.rho(p_A,T_BED1max);
    rho_Amax=DryAir.rho(p_A,T_BED1min);

    Vdotmin=mean([VdotAir15min,VdotAir20min],2).*rho_Anorm./rho_Amax;
    Vdotmax=mean([VdotAir15max,VdotAir20max],2).*rho_Anorm./rho_Amin;
    
    w_mfMin=w_mf(d_p,rho_p,p_A,T_BED1max);
    w_mfMax=w_mf(d_p,rho_p,p_A,T_BED1min);
    
    Vdot_mfMin=w_mfMin.*l.*w;
    Vdot_mfMax=w_mfMax.*l.*w;

    tab.FGlower=Vdotmin./Vdot_mfMax;
    tab.FGupper=Vdotmax./Vdot_mfMin;


    %% Heat transfer coefficients
    T_surfMin=mean([T_PROBE1min,T_PROBE2min],2);
    T_surfMax=mean([T_PROBE1max,T_PROBE2max],2);

    tab.alpha_grossLower=Umin.*Imin./(A.*(T_surfMax-T_BED1min));
    tab.alpha_grossUpper=Umax.*Imax./(A.*(T_surfMin-T_BED1max));


    mDotAmin=mean([VdotAir15min,VdotAir20min],2).*rho_Anorm;
    mDotAmax=mean([VdotAir15max,VdotAir20max],2).*rho_Anorm;

    QdotLossMin=mDotAmin.*(DryAir.h(T_LausMin)-DryAir.h(T_LeinMax));
    QdotLossMax=mDotAmax.*(DryAir.h(T_LausMax)-DryAir.h(T_LeinMin));

    tab.alpha_netLower=(Umin.*Imin-QdotLossMax)./(A.*(T_surfMax-T_BED1min));
    tab.alpha_netUpper=(Umax.*Imax-QdotLossMin)./(A.*(T_surfMin-T_BED1max));
end


function [Tmin,Tmax]=TRange(T)
    delta=0.6+0.0017*abs(T-273.15);
    
    Tmin=T-delta;
    Tmax=T+delta;
end


function [Vdotmin,Vdotmax]=VdotRange(Vdot,VdotNorm)
    VdotNorm=VdotNorm.*60^2;
    
    delta=NaN(size(Vdot));
    delta(VdotNorm<10)=1e-2*10./VdotNorm(VdotNorm<10);
    delta(VdotNorm>=10)=1e-2;
    delta=delta.*Vdot;
    
    Vdotmin=Vdot-delta;
    Vdotmax=Vdot+delta;
end


function [Umin,Umax]=URange(U)
    delta=0.05e-2.*U+20e-3;
    
    Umin=U-delta;
    Umax=U+delta;
end


function [Imin,Imax]=IRange(I)
    delta=0.2e-2.*I+20e-3;
    
    Imin=I-delta;
    Imax=I+delta;
end






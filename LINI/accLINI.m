%% Accuracy of LINI test rig results
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
%This function calculates the measurement uncertainties of all LINI test 
%rig results.
%This function is only called during the general analysis of the results in
%the script "Analyze_LINI.m".
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir
%   - w_mf.m


function tab=accLINI(tab)
    %% Constants
    l_in=0.15;
    l_out=l_in;
    l_1=0.293;
    l_2=l_1;

    l_heated=0.25;
    d=25e-3;
    A=d*pi*l_heated;

    A_FB=0.9*0.194;
    x_heated=l_1/(l_in+l_1+l_2+l_out);

    DeltaH_eps=50e-3;

    d_p=146e-6;
    rho_p=2650;
    rho_Anorm=1.293;
    g=9.81;
    p_amb=1013.25e2;
    
    
    %% Calculate ranges
    [p_eps1min,p_eps1max]=pRange(tab.p_eps1,0,10e2);
    [p_eps2min,p_eps2max]=pRange(tab.p_eps2,0,10e2);
    [dpFloormin,dpFloormax]=pRange(tab.dpFloor,0,25e2);
    [pAinmin,pAinmax]=pRange(tab.pAin,0,250e2);    
    
    [TAinMin,TAinMax]=TRange(tab.TAin);
    [TAoutMin,TAoutMax]=TRange(tab.TAout);
    [T_bedMin,T_bedMax]=TRange(tab.T_bed);
    [T_surfMin,T_surfMax]=TRange(tab.T_surf);

    [VdotRotaMin,VdotRotaMax]=VdotRange(tab.Vdot,tab.VdotRota);

    [Umin,Umax]=URange(tab.U);
    [Imin,Imax]=IRange(tab.I);
    [phimin,phimax]=phiRange(tab.phi);


    %% Degree of Fluidization
    %Calculate bed porosities
    eps1min=1-p_eps1max./(rho_p*g*DeltaH_eps);
    eps1max=1-p_eps1min./(rho_p*g*DeltaH_eps);
    
    eps2min=1-p_eps2max./(rho_p*g*DeltaH_eps);
    eps2max=1-p_eps2min./(rho_p*g*DeltaH_eps);
    
    epsmin=mean([eps1min,eps2min],2);
    epsmax=mean([eps1max,eps2max],2);
    
    
    p_Amin=p_amb+pAinmin-dpFloormax+rho_p*g*35e-3*(1-epsmax);
    p_Amax=p_amb+pAinmax-dpFloormin+rho_p*g*35e-3*(1-epsmin);
    
    T_Amin=mean([TAinMin,TAoutMin],2);
    T_Amax=mean([TAinMax,TAoutMax],2);
    
    rho_Amin=DryAir.rho(p_Amin,T_Amax);
    rho_Amax=DryAir.rho(p_Amax,T_Amin);

    VdotMin=VdotRotaMin.*rho_Anorm./rho_Amax;
    VdotMax=VdotRotaMax.*rho_Anorm./rho_Amin;
    
    w_mfMin=w_mf(d_p,rho_p,p_Amin,T_Amax);
    w_mfMax=w_mf(d_p,rho_p,p_Amax,T_Amin);
    
    Vdot_mfMin=w_mfMin.*A_FB;
    Vdot_mfMax=w_mfMax.*A_FB;

    tab.FGlower=VdotMin./Vdot_mfMax;
    tab.FGupper=VdotMax./Vdot_mfMin;


    %% Heat transfer coefficients
    tab.alpha_grossLower=Umin.*Imin.*cos(phimax)./(A.*(T_surfMax-T_bedMin));
    tab.alpha_grossUpper=Umax.*Imax.*cos(phimin)./(A.*(T_surfMin-T_bedMax));

    
    mDotAmin=VdotRotaMin.*rho_Anorm.*x_heated;
    mDotAmax=VdotRotaMax.*rho_Anorm.*x_heated;

    QdotLossMin=mDotAmin.*(DryAir.h(TAoutMin)-DryAir.h(TAinMax));
    QdotLossMax=mDotAmax.*(DryAir.h(TAoutMax)-DryAir.h(TAinMin));

    tab.alpha_netLower=(Umin.*Imin.*cos(phimax)-QdotLossMax)./(A.*(T_surfMax-T_bedMin));
    tab.alpha_netUpper=(Umax.*Imax.*cos(phimin)-QdotLossMin)./(A.*(T_surfMin-T_bedMax));
end


function [Tmin,Tmax]=TRange(T)
    persistent A B R0 offset gain
    if isempty(A)
        A=3.9083e-3;
        B=-5.775e-7;
        R0=100;
        
        offset=0.0015e-2*(390-0.5)/(390.48-18.52)*(850--200);
        gain=0.0059e-2;
    end
    T=T-273.15;
    
    delta=0.15+0.002*abs(T);
    Tmin=T-delta;
    Tmax=T+delta;
    
    Rmin=Rfx(Tmin)*(1-gain);
    Rmax=Rfx(Tmax)*(1+gain);
    
    for i=1:length(T)
        Tmin(i)=fzero(@(T) Rfx(T)-Rmin(i),[Tmin(i)-10,Tmin(i)+10]);
        Tmax(i)=fzero(@(T) Rfx(T)-Rmax(i),[Tmax(i)-10,Tmax(i)+10]);
    end
    
    Tmin=Tmin-offset+273.15;
    Tmax=Tmax+offset+273.15;
    
    
    function R=Rfx(T)
        R=R0*(1+A*T+B*T.^2);
    end
end


function [Vdotmin,Vdotmax]=VdotRange(Vdot,VdotNorm)
    VdotNorm=VdotNorm.*60^2;
    
    delta=NaN(size(Vdot));
    delta(VdotNorm<22)=1e-2*22./VdotNorm(VdotNorm<22);
    delta(VdotNorm>=22)=1e-2;
    delta=delta.*Vdot;
    
    Vdotmin=Vdot-delta;
    Vdotmax=Vdot+delta;
end


function [pmin,pmax]=pRange(p,MRB,MRE)
    range=MRE-MRB;
    if range<=10e2
        bias=1e-2*range;
    elseif range<=100e2
        bias=0.8e-2*range;
    else
        bias=0.5e-2*range;
    end
    
    
    pmin=X20AI4632min(p-bias,MRB,MRE);
    pmax=X20AI4632max(p+bias,MRB,MRE);
end


function [Umin,Umax]=URange(U)
    delta=0.65e-2.*U;
    
    Umin=U-delta;
    Umax=U+delta;
end


function [Imin,Imax]=IRange(I)
    delta=1.707e-2.*I;
    
    Imin=I-delta;
    Imax=I+delta;
end


function [phimin,phimax]=phiRange(phi)
    phi=rad2deg(phi);
    
    gain=0.5e-2;
    bias=1;
    
    phimin=(phi-bias).*(1-gain);
    phimax=(phi+bias).*(1+gain);
    
    phimin(phimin<-90)=-90;
    phimax(phimax>90)=90;
    
    phimin=deg2rad(phimin);
    phimax=deg2rad(phimax);
end


function xmin=X20AI4632min(x,MRB,MRE)
    persistent gain bias
    if isempty(gain)
        gain=0.08e-2;
        bias=0.02e-2;
    end
    
    xmin=((x./(MRE-MRB)*16+4)*(1-gain)-bias*20-4)/16*(MRE-MRB);
end


function xmax=X20AI4632max(x,MRB,MRE)
    persistent gain bias
    if isempty(gain)
        gain=0.08e-2;
        bias=0.02e-2;
    end
    
    xmax=((x./(MRE-MRB)*16+4)*(1+gain)+bias*20-4)/16*(MRE-MRB);
end






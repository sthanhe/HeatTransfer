%% Accuracy of TRINI test rig results
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
%This function calculates the measurement uncertainties of all TRINI test 
%rig results.
%This function is only called during the general analysis of the results in
%the script "Analyze_TRINI.m".
%
%Required products:
%   - MATLAB, version 9.10
%Necessary files, classes and functions:
%   - @DryAir
%   - w_mf.m


function tab=accTRINI(tab,plain)
    %% Constants
    l_heated=0.216;
    d=0.025;
    A=d*pi*l_heated;

    A_FBin=0.2*0.203;
    A_FBout=A_FBin;
    A_FBmain=0.2*0.597;

    d_p=146e-6;
    rho_p=2650;

    eta_A=18.107811e-6;
    R_A=287.0533;

    p_amb=101325;

    alpha=10e-12;
    beta=83e-12;
    s=20e-3;
    
    
    %% Calculate ranges
    [p1min,p1max]=pRange(tab.p1,0,50e3);
    [p2min,p2max]=pRange(tab.p2,0,5e3);
    [p3min,p3max]=pRange(tab.p3,0,5e3);
    [p4min,p4max]=pRange(tab.p4,0,25e3);
    [p5min,p5max]=pRange(tab.p5,0,25e3);
    [p6min,p6max]=pRange(tab.p6,0,25e3);
    [p7min,p7max]=pRange(tab.p7,0,25e3);
    [p8min,p8max]=pRange(tab.p8,0,25e3);
    
    [T2min,T2max]=TRange(tab.T2);
    [T3min,T3max]=TRange(tab.T3);
    [T4min,T4max]=TRange(tab.T4);
    [TAoutMin,TAoutMax]=TRange(tab.TAout);
    [T_surf1min,T_surf1max]=TRange(tab.T_surf1);
    [T_surf2min,T_surf2max]=TRange(tab.T_surf2);
    [T_bed11min,T_bed11max]=TRange(tab.T_bed11);
    [T_bed12min,T_bed12max]=TRange(tab.T_bed12);
    [T_bed21min,T_bed21max]=TRange(tab.T_bed21);
    [T_bed22min,T_bed22max]=TRange(tab.T_bed22);
    
    if ~plain
        [T_bed13min,T_bed13max]=TRange(tab.T_bed13);
        [T_bed14min,T_bed14max]=TRange(tab.T_bed14);
        [T_bed23min,T_bed23max]=TRange(tab.T_bed23);
        [T_bed24min,T_bed24max]=TRange(tab.T_bed24);
    end
    
    [U1min,U1max]=URange(tab.U1);
    [U2min,U2max]=URange(tab.U2);
    
    [I1min,I1max]=IRange(tab.I1);
    [I2min,I2max]=IRange(tab.I2);
    
    [phi1min,phi1max]=phiRange(tab.phi1);
    [phi2min,phi2max]=phiRange(tab.phi2);
    
    [mDotAmin,mDotAmax]=mDotRange(tab.mDotA);


    %% Degree of Fluidization
    pInMin=p1min-p2max+p_amb;
    pInMax=p1max-p2min+p_amb;
    
    pMainMin=p1min-p3max+p_amb;
    pMainMax=p1max-p3min+p_amb;
    
    pOutMin=p1min-p4max+p_amb;
    pOutMax=p1max-p4min+p_amb;
    

    rhoInMin=pInMin./(R_A.*T2max);
    rhoInMax=pInMax./(R_A.*T2min);
    
    rhoMainMin=pMainMin./(R_A.*T3max);
    rhoMainMax=pMainMax./(R_A.*T3min);
    
    rhoOutMin=pOutMin./(R_A.*T4max);
    rhoOutMax=pOutMax./(R_A.*T4min);
    

    lambdaInMin=eta_A.*beta.*A_FBin./(rhoInMax.*alpha);
    lambdaInMax=eta_A.*beta.*A_FBin./(rhoInMin.*alpha);
    
    lambdaMainMin=eta_A.*beta.*A_FBmain./(rhoMainMax.*alpha);
    lambdaMainMax=eta_A.*beta.*A_FBmain./(rhoMainMin.*alpha);
    
    lambdaOutMin=eta_A.*beta.*A_FBout./(rhoOutMax.*alpha);
    lambdaOutMax=eta_A.*beta.*A_FBout./(rhoOutMin.*alpha);
    

    VdotIn_estMin=-lambdaInMin./2+sqrt((lambdaInMin./2).^2+beta.*A_FBin.^2.*p5min./(rhoInMax.*s));
    VdotIn_estMax=-lambdaInMax./2+sqrt((lambdaInMax./2).^2+beta.*A_FBin.^2.*p5max./(rhoInMin.*s));
    
    VdotMain_estMin=-lambdaMainMin./2+sqrt((lambdaMainMin./2).^2+beta.*A_FBmain.^2.*(p6min+p7min)./(2*rhoMainMax.*s));
    VdotMain_estMax=-lambdaMainMax./2+sqrt((lambdaMainMax./2).^2+beta.*A_FBmain.^2.*(p6max+p7max)./(2*rhoMainMin.*s));
    
    VdotOut_estMin=-lambdaOutMin./2+sqrt((lambdaOutMin./2).^2+beta.*A_FBout.^2.*p8min./(rhoOutMax.*s));
    VdotOut_estMax=-lambdaOutMax./2+sqrt((lambdaOutMax./2).^2+beta.*A_FBout.^2.*p8max./(rhoOutMin.*s));
    

    mDotAin_estMin=VdotIn_estMin.*rhoInMax;
    mDotAin_estMax=VdotIn_estMax.*rhoInMin;
    
    mDotAmain_estMin=VdotMain_estMin.*rhoMainMax;
    mDotAmain_estMax=VdotMain_estMax.*rhoMainMin;
    
    mDotAout_estMin=VdotOut_estMin.*rhoOutMax;
    mDotAout_estMax=VdotOut_estMax.*rhoOutMin;
    

    mDotA_estMin=mDotAin_estMin+mDotAmain_estMin+mDotAout_estMin;
    mDotA_estMax=mDotAin_estMax+mDotAmain_estMax+mDotAout_estMax;

    
%     mDotAinMin=mDotAin_estMin./mDotA_estMin.*mDotAmin;
%     mDotAinMax=mDotAin_estMax./mDotA_estMax.*mDotAmax;
    
    mDotAmainMin=mDotAmain_estMin./mDotA_estMin.*mDotAmin;
    mDotAmainMax=mDotAmain_estMax./mDotA_estMax.*mDotAmax;
    
%     mDotAoutMin=mDotAout_estMin./mDotA_estMin.*mDotAmin;
%     mDotAoutMax=mDotAout_estMax./mDotA_estMax.*mDotAmax;

    pAinMin=pMainMin-mean([p6max,p7max],2);
    pAinMax=pMainMax-mean([p6min,p7min],2);

    TAinMin=T3min;
    TAinMax=T3max;
    
    T_Amin=mean([TAinMin,TAoutMin],2);
    T_Amax=mean([TAinMax,TAoutMax],2);
    
    rho_Amin=DryAir.rho(pAinMin,T_Amax);
    rho_Amax=DryAir.rho(pAinMax,T_Amin);
    
    w_mfMin=w_mf(d_p,rho_p,pAinMin,T_Amax);
    w_mfMax=w_mf(d_p,rho_p,pAinMax,T_Amin);
    
    Vdot_mfMainMin=w_mfMin.*A_FBmain;
    Vdot_mfMainMax=w_mfMax.*A_FBmain;
    
    tab.FGlower=mDotAmainMin./(rho_Amax.*Vdot_mfMainMax);
    tab.FGupper=mDotAmainMax./(rho_Amin.*Vdot_mfMainMin);


    %% Heat transfer coefficients
    %Gross
    if plain
        T_bed1min=mean([T_bed11min,T_bed12min],2);
        T_bed1max=mean([T_bed11max,T_bed12max],2);
        
        T_bed2min=mean([T_bed21min,T_bed22min],2);
        T_bed2max=mean([T_bed21max,T_bed22max],2);
    else
        T_bed1min=mean([T_bed11min,T_bed12min,T_bed13min,T_bed14min],2);
        T_bed1max=mean([T_bed11max,T_bed12max,T_bed13max,T_bed14max],2);
        
        T_bed2min=mean([T_bed21min,T_bed22min,T_bed23min,T_bed24min],2);
        T_bed2max=mean([T_bed21max,T_bed22max,T_bed23max,T_bed24max],2);
    end
    

    alpha_gross1Min=U1min.*I1min.*cos(phi1max)./(A.*(T_surf1max-T_bed1min));
    alpha_gross1Max=U1max.*I1max.*cos(phi1min)./(A.*(T_surf1min-T_bed1max));
    
    alpha_gross2Min=U2min.*I2min.*cos(phi2max)./(A.*(T_surf2max-T_bed2min));
    alpha_gross2Max=U2max.*I2max.*cos(phi2min)./(A.*(T_surf2min-T_bed2max));
    
    tab.alpha_grossMeanLower=mean([alpha_gross1Min,alpha_gross2Min],2);
    tab.alpha_grossMeanUpper=mean([alpha_gross1Max,alpha_gross2Max],2);


    %Net
    QdotLossMin=mDotAmainMin.*(DryAir.h(TAoutMin)-DryAir.h(TAinMax));
    QdotLossMax=mDotAmainMax.*(DryAir.h(TAoutMax)-DryAir.h(TAinMin));
    
    alpha_net1min=(U1min.*I1min.*cos(phi1max)-QdotLossMax./2)./A./(T_surf1max-T_bed1min);
    alpha_net1max=(U1max.*I1max.*cos(phi1min)-QdotLossMin./2)./A./(T_surf1min-T_bed1max);
    
    alpha_net2min=(U2min.*I2min.*cos(phi2max)-QdotLossMax./2)./A./(T_surf2max-T_bed2min);
    alpha_net2max=(U2max.*I2max.*cos(phi2min)-QdotLossMin./2)./A./(T_surf2min-T_bed2max);
    
    tab.alpha_netMeanLower=mean([alpha_net1min,alpha_net2min],2);
    tab.alpha_netMeanUpper=mean([alpha_net1max,alpha_net2max],2);
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


function [mDotmin,mDotmax]=mDotRange(mDot)
    mDot100=2030/60^2;
    lim=0.2*mDot100;
    
    MRB=0;
    MRE=0.5;
    biasTransformer=10e-3/16*(MRE-MRB);
    
    delta=NaN(size(mDot));
    delta(mDot<lim)=0.8e-2*mDot100;
    delta(mDot>=lim)=4e-2*mDot(mDot>=lim);
    
    mDotmin=mDot-delta-biasTransformer;
    mDotmax=mDot+delta+biasTransformer;
    
    mDotmin=X20AI4632min(mDotmin,MRB,MRE);
    mDotmax=X20AI4632max(mDotmax,MRB,MRE);
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






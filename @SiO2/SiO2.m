%% Property functions of silicon dioxide (SiO2)
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
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5500329
%
%
%
%This class describes the thermo-physical properties of silicon dioxide
%(SiO2) according to:
%
%Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. 
%Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951
%https://webbook.nist.gov/cgi/cbook.cgi?ID=C14808607&Type=JANAFS&Table=on
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - createFits.m
%   - fits.mat


classdef SiO2
    %All parameters and results in SI base units
    
    %%
    properties(Constant)
        M=60.0843e-3;   %Molar Mass
    end
    
    properties(Constant, Access=private)
        A=-6.076591;
        B=251.6755;
        C=-324.7964;
        D=168.5604;
        E=0.002548;
        F=-917.6893;

        A_beta=58.75340;
        B_beta=10.27925;
        C_beta=-0.131384;
        D_beta=0.025210;
        E_beta=0.025601;
        F_beta=-929.3292;
    end
    
    
    %% Property Functions
    methods(Static)
        function c_p=c_p(T)
            %Specific isobaric heat capacity
            [alpha,beta]=SiO2.getPhases(T);
            T=T./1000;
            
            c_p=NaN(size(T));
            c_p(alpha)=c_pfx(T(alpha),SiO2.A,SiO2.B,SiO2.C,SiO2.D,SiO2.E);
            c_p(beta)=c_pfx(T(beta),SiO2.A_beta,SiO2.B_beta,SiO2.C_beta,SiO2.D_beta,SiO2.E_beta);
            
            
            function c_p=c_pfx(T,A,B,C,D,E)
                c_p=(A+B*T+C*T.^2+D*T.^3+E./T.^2)./SiO2.M;
            end
        end
        
        
        function h=h(T)
            %Specific enthalpy
            %h(298.15)=0
            persistent H
            if isempty(H)
                H=-910.8568;
            end
            
            [alpha,beta]=SiO2.getPhases(T);
            T=T./1000;
            
            h=NaN(size(T));
            h(alpha)=hfx(T(alpha),SiO2.A,SiO2.B,SiO2.C,SiO2.D,SiO2.E,SiO2.F);
            h(beta)=hfx(T(beta),SiO2.A_beta,SiO2.B_beta,SiO2.C_beta,SiO2.D_beta,SiO2.E_beta,SiO2.F_beta);
            
            
            function h=hfx(T,A,B,C,D,E,F)
                h=(A*T+B*T.^2./2+C*T.^3./3+D*T.^4./4-E./T+F-H)./SiO2.M.*1000;
            end
        end
        
        
        function s=s(T)
            %Specific entropy
            persistent G G_beta
            if isempty(G)
                G=-27.96962;
                G_beta=105.8092;
            end
            
            [alpha,beta]=SiO2.getPhases(T);
            T=T./1000;
            
            s=NaN(size(T));
            s(alpha)=sfx(T(alpha),SiO2.A,SiO2.B,SiO2.C,SiO2.D,SiO2.E,G);
            s(beta)=sfx(T(beta),SiO2.A_beta,SiO2.B_beta,SiO2.C_beta,SiO2.D_beta,SiO2.E_beta,G_beta);
            
            
            function s=sfx(T,A,B,C,D,E,G)
                s=(A*log(T)+B*T+C*T.^2./2+D*T.^3./3-E./(2*T.^2)+G)./SiO2.M;
            end
        end
    end
    
    
    %% Backward Equations
    methods(Static)
        function T=T_h(h)
            %Backwards-equation for temperature as function of specific
            %enthalpy
            persistent T_halpha T_hbeta hmin hbound hmax
            if isempty(T_halpha)
                vars=load('SiO2\fits.mat','T_halpha','T_hbeta');
                T_halpha=vars.T_halpha;
                T_hbeta=vars.T_hbeta;
                
                hmin=SiO2.h(298);
                hbound=SiO2.h(847);
                hmax=SiO2.h(1996);
            end
            
            alpha=hmin<=h & h<hbound;
            beta=hbound<=h & h<=hmax;
            
            T=NaN(size(h));
            T(alpha)=T_halpha(h(alpha));
            T(beta)=T_hbeta(h(beta));
        end
        
        
        function createConstants()
            %Creates the curve fittings for lookup
            clear('SiO2');
            
            Talpha=298:0.01:847-0.01;
            halpha=SiO2.h(Talpha);
            Tbeta=847:0.01:1996;
            hbeta=SiO2.h(Tbeta);
            
            
            fitresult=SiO2.createFits(halpha,Talpha,hbeta,Tbeta);
            
            T_halpha=fitresult{1};
            T_hbeta=fitresult{2};
            
            save('fits.mat','T_halpha','T_hbeta');
        end
    end
    
    
    %% Internal Functions
    methods(Static, Access=private)
        [fitresult,gof]=createFits(halpha,Talpha,hbeta,Tbeta)
        
        
        function [alpha,beta]=getPhases(T)
            alpha=298<=T & T<847;
            beta=847<=T & T<=1996;
        end
    end
end





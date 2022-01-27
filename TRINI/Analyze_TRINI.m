%% Analyze TRINI test rig results
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
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5500329
%
%
%
%This script analyzes the TRINI test rig data and creates all published
%figures.
%
%Requires all TRINI data files in the same folder and all auxiliary
%functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
%   - accTRINI.m
%   - w_mf.m


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


dirCont=dir();


%% Plain tubes, without baffle
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_PlainWithout_'));


%Initialize table
varnames={'Dataset','I1','I2','T2','T3',...
            'T4','TAout','T_bed11','T_bed12','T_bed13',...
            'T_bed14','T_bed21','T_bed22','T_bed23','T_bed24',...
            'T_surf1','T_surf2','U1','U2','mDotA','p1',...
            'p2','p3','p4','p5','p6',...
            'p7','p8','phi1','phi2',...
            'FG','FGlower','FGupper',...
            'alpha_gross1','alpha_gross2','alpha_net1','alpha_net2',...
            'alpha_grossMean','alpha_grossMeanLower','alpha_grossMeanUpper',...
            'alpha_netMean','alpha_netMeanLower','alpha_netMeanUpper'};
tabPlainWithout=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tabPlainWithout.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tabPlainWithout{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:13}=tabloc{:,4:13}+273.15;
    tabloc{:,25:26}=deg2rad(tabloc{:,25:26});
    
    tabPlainWithout{i,[2:9,12:13,16:30]}=mean(tabloc{:,2:end});
    tabPlainWithout{i,[10,11,14,15]}=NaN;
end


%Calculate degree of fluidization
pIn=tabPlainWithout.p1-tabPlainWithout.p2+p_amb;
pMain=tabPlainWithout.p1-tabPlainWithout.p3+p_amb;
pOut=tabPlainWithout.p1-tabPlainWithout.p4+p_amb;

rhoIn=pIn./(R_A.*(tabPlainWithout.T2));
rhoMain=pMain./(R_A.*(tabPlainWithout.T3));
rhoOut=pOut./(R_A.*(tabPlainWithout.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tabPlainWithout.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tabPlainWithout.p6+tabPlainWithout.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tabPlainWithout.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tabPlainWithout.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

p_AplainWithout=pMain-mean([tabPlainWithout.p6,tabPlainWithout.p7],2);
TAin=tabPlainWithout.T3;
T_AplainWithout=mean([TAin,tabPlainWithout.TAout],2);
rho_A=DryAir.rho(p_AplainWithout,T_AplainWithout);
wmfplainWithout=w_mf(d_p,rho_p,p_AplainWithout,T_AplainWithout);
Vdot_mfMain=wmfplainWithout.*A_FBmain;
tabPlainWithout.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FGwithout=tabPlainWithout.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tabPlainWithout.T_bed11,tabPlainWithout.T_bed12],2);
T_bed2=mean([tabPlainWithout.T_bed21,tabPlainWithout.T_bed22],2);

tabPlainWithout.alpha_gross1=tabPlainWithout.U1.*tabPlainWithout.I1.*cos(tabPlainWithout.phi1)./A./(tabPlainWithout.T_surf1-T_bed1);
tabPlainWithout.alpha_gross2=tabPlainWithout.U2.*tabPlainWithout.I2.*cos(tabPlainWithout.phi2)./A./(tabPlainWithout.T_surf2-T_bed2);
tabPlainWithout.alpha_grossMean=mean([tabPlainWithout.alpha_gross1,tabPlainWithout.alpha_gross2],2);
alpha_grossMeanWithout=tabPlainWithout.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tabPlainWithout.TAout)-DryAir.h(TAin));
tabPlainWithout.alpha_net1=(tabPlainWithout.U1.*tabPlainWithout.I1.*cos(tabPlainWithout.phi1)-QdotLoss./2)./A./(tabPlainWithout.T_surf1-T_bed1);
tabPlainWithout.alpha_net2=(tabPlainWithout.U2.*tabPlainWithout.I2.*cos(tabPlainWithout.phi2)-QdotLoss./2)./A./(tabPlainWithout.T_surf2-T_bed2);
tabPlainWithout.alpha_netMean=mean([tabPlainWithout.alpha_net1,tabPlainWithout.alpha_net2],2);
alpha_netMeanWithout=tabPlainWithout.alpha_netMean;


%Calculate accuracy
tabPlainWithout=accTRINI(tabPlainWithout,true);


%Create figure
fig=figure(12);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6.5];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=tabPlainWithout.FG;
FGbounds=[1.5,2,2.5,3.3,3.6,4,4.5,5.2,5.7,6.5];
[gofGross1,gofNet1]=plotGrossNet(ax,FG,FGbounds,tabPlainWithout.alpha_gross1,tabPlainWithout.alpha_net1,xfit,[1,2]);
[gofGross2,gofNet2]=plotGrossNet(ax,FG,FGbounds,tabPlainWithout.alpha_gross2,tabPlainWithout.alpha_net2,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['HTCs TRINI Rig, Plain without Baffle, w_{mf}=',num2str(round(mean(wmfplainWithout*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{['Gross1, R^2 = ',num2str(round(gofGross1.rsquare,3))],...
                    ['Net1, R^2 = ',num2str(round(gofNet1.rsquare,3))],...
                    ['Gross2, R^2 = ',num2str(round(gofGross2.rsquare,3))],...
                    ['Net2, R^2 = ',num2str(round(gofNet2.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';

ax2=axes(fig,'Units','centimeters','Position',[2,5,3,3]);
imshow('withoutBaffle.tiff','Parent',ax2,'Reduce',false);

saveas(fig,'Figure12.tiff');

save('resultsTRINI.mat','tabPlainWithout','p_AplainWithout','T_AplainWithout','wmfplainWithout');


%% Plain tubes, with baffle
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_PlainWith_'));


%Initialize table
tabPlainWith=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tabPlainWith.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tabPlainWith{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:13}=tabloc{:,4:13}+273.15;
    tabloc{:,25:26}=deg2rad(tabloc{:,25:26});
    
    tabPlainWith{i,[2:9,12:13,16:30]}=mean(tabloc{:,2:end});
    tabPlainWith{i,[10,11,14,15]}=NaN;
end


%Calculate degree of fluidization
pIn=tabPlainWith.p1-tabPlainWith.p2+p_amb;
pMain=tabPlainWith.p1-tabPlainWith.p3+p_amb;
pOut=tabPlainWith.p1-tabPlainWith.p4+p_amb;

rhoIn=pIn./(R_A.*(tabPlainWith.T2));
rhoMain=pMain./(R_A.*(tabPlainWith.T3));
rhoOut=pOut./(R_A.*(tabPlainWith.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tabPlainWith.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tabPlainWith.p6+tabPlainWith.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tabPlainWith.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tabPlainWith.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

p_AplainWith=pMain-mean([tabPlainWith.p6,tabPlainWith.p7],2);
TAin=tabPlainWith.T3;
T_AplainWith=mean([TAin,tabPlainWith.TAout],2);
rho_A=DryAir.rho(p_AplainWith,T_AplainWith);
wmfplainWith=w_mf(d_p,rho_p,p_AplainWith,T_AplainWith);
Vdot_mfMain=wmfplainWith.*A_FBmain;
tabPlainWith.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FGwith=tabPlainWith.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tabPlainWith.T_bed11,tabPlainWith.T_bed12],2);
T_bed2=mean([tabPlainWith.T_bed21,tabPlainWith.T_bed22],2);

tabPlainWith.alpha_gross1=tabPlainWith.U1.*tabPlainWith.I1.*cos(tabPlainWith.phi1)./A./(tabPlainWith.T_surf1-T_bed1);
tabPlainWith.alpha_gross2=tabPlainWith.U2.*tabPlainWith.I2.*cos(tabPlainWith.phi2)./A./(tabPlainWith.T_surf2-T_bed2);
tabPlainWith.alpha_grossMean=mean([tabPlainWith.alpha_gross1,tabPlainWith.alpha_gross2],2);
alpha_grossMeanWith=tabPlainWith.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tabPlainWith.TAout)-DryAir.h(TAin));
tabPlainWith.alpha_net1=(tabPlainWith.U1.*tabPlainWith.I1.*cos(tabPlainWith.phi1)-QdotLoss./2)./A./(tabPlainWith.T_surf1-T_bed1);
tabPlainWith.alpha_net2=(tabPlainWith.U2.*tabPlainWith.I2.*cos(tabPlainWith.phi2)-QdotLoss./2)./A./(tabPlainWith.T_surf2-T_bed2);
tabPlainWith.alpha_netMean=mean([tabPlainWith.alpha_net1,tabPlainWith.alpha_net2],2);
alpha_netMeanWith=tabPlainWith.alpha_netMean;


%Calculate accuracy
tabPlainWith=accTRINI(tabPlainWith,true);


%Create figure
fig=figure(13);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6.5];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=tabPlainWith.FG;
FGbounds=[1.5,2.2,2.7,3.5,4.5,6];
[gofGross1,gofNet1]=plotGrossNet(ax,FG,FGbounds,tabPlainWith.alpha_gross1,tabPlainWith.alpha_net1,xfit,[1,2]);
[gofGross2,gofNet2]=plotGrossNet(ax,FG,FGbounds,tabPlainWith.alpha_gross2,tabPlainWith.alpha_net2,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['HTCs TRINI Rig, Plain with Baffle, w_{mf}=',num2str(round(mean(wmfplainWith*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{['Gross1, R^2 = ',num2str(round(gofGross1.rsquare,3))],...
                    ['Net1, R^2 = ',num2str(round(gofNet1.rsquare,3))],...
                    ['Gross2, R^2 = ',num2str(round(gofGross2.rsquare,3))],...
                    ['Net2, R^2 = ',num2str(round(gofNet2.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';

ax2=axes(fig,'Units','centimeters','Position',[2,5,3,3]);
imshow('withBaffle.tiff','Parent',ax2);

saveas(fig,'Figure13.tiff');

save('resultsTRINI.mat','tabPlainWith','p_AplainWith','T_AplainWith','wmfplainWith','-append');


%% Plain tubes, comparison with / without baffle
%Create figure
fig=figure(14);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6.5];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=[tabPlainWith.FG;tabPlainWith.FG];
alpha_gross=[tabPlainWith.alpha_gross1;tabPlainWith.alpha_gross2];
alpha_net=[tabPlainWith.alpha_net1;tabPlainWith.alpha_net2];
FGbounds=[1.5,2.2,2.7,3.5,4.5,6];
[gofGrossWith,gofNetWith]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,[1,2]);


FG=[tabPlainWithout.FG;tabPlainWithout.FG];
alpha_gross=[tabPlainWithout.alpha_gross1;tabPlainWithout.alpha_gross2];
alpha_net=[tabPlainWithout.alpha_net1;tabPlainWithout.alpha_net2];
FGbounds=[1.5,2,2.5,3.3,3.6,4,4.5,5.2,5.7,6.5];
[gofGrossWithout,gofNetWithout]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'HTCs TRINI Rig, Plain with and without Baffle');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,legItems,{['Gross With, R^2 = ',num2str(round(gofGrossWith.rsquare,3))],...
                    ['Net With, R^2 = ',num2str(round(gofNetWith.rsquare,3))],...
                    ['Gross Without, R^2 = ',num2str(round(gofGrossWithout.rsquare,3))],...
                    ['Net Without, R^2 = ',num2str(round(gofNetWithout.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure14.tiff');


%% Finned tubes, fin pitch 9 mm, fin thickness 2 mm (9/2)
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_Finned92_'));


%Initialize table
tab92=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab92.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab92{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:17}=tabloc{:,4:17}+273.15;
    tabloc{:,29:30}=deg2rad(tabloc{:,29:30});
    
    tab92{i,2:30}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab92.p1-tab92.p2+p_amb;
pMain=tab92.p1-tab92.p3+p_amb;
pOut=tab92.p1-tab92.p4+p_amb;

rhoIn=pIn./(R_A.*(tab92.T2));
rhoMain=pMain./(R_A.*(tab92.T3));
rhoOut=pOut./(R_A.*(tab92.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab92.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab92.p6+tab92.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab92.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab92.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

p_A92=pMain-mean([tab92.p6,tab92.p7],2);
TAin=tab92.T3;
T_A92=mean([TAin,tab92.TAout],2);
rho_A=DryAir.rho(p_A92,T_A92);
wmf92=w_mf(d_p,rho_p,p_A92,T_A92);
Vdot_mfMain=wmf92.*A_FBmain;
tab92.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FG92=tab92.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tab92.T_bed11,tab92.T_bed12,tab92.T_bed13,tab92.T_bed14],2);
T_bed2=mean([tab92.T_bed21,tab92.T_bed22,tab92.T_bed23,tab92.T_bed24],2);

tab92.alpha_gross1=tab92.U1.*tab92.I1.*cos(tab92.phi1)./A./(tab92.T_surf1-T_bed1);
tab92.alpha_gross2=tab92.U2.*tab92.I2.*cos(tab92.phi2)./A./(tab92.T_surf2-T_bed2);
tab92.alpha_grossMean=mean([tab92.alpha_gross1,tab92.alpha_gross2],2);
alpha_grossMean92=tab92.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab92.TAout)-DryAir.h(TAin));
tab92.alpha_net1=(tab92.U1.*tab92.I1.*cos(tab92.phi1)-QdotLoss./2)./A./(tab92.T_surf1-T_bed1);
tab92.alpha_net2=(tab92.U2.*tab92.I2.*cos(tab92.phi2)-QdotLoss./2)./A./(tab92.T_surf2-T_bed2);
tab92.alpha_netMean=mean([tab92.alpha_net1,tab92.alpha_net2],2);
alpha_netMean92=tab92.alpha_netMean;


%Calculate accuracy
tab92=accTRINI(tab92,false);


%Create figure
fig=figure(15);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[2,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=tab92.FG;
FGbounds=[2,2.5,3,3.7,4.5,6];
[gofGross1,gofNet1]=plotGrossNet(ax,FG,FGbounds,tab92.alpha_gross1,tab92.alpha_net1,xfit,[1,2]);
[gofGross2,gofNet2]=plotGrossNet(ax,FG,FGbounds,tab92.alpha_gross2,tab92.alpha_net2,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['Virtual HTCs TRINI Rig, Finned 9/2, w_{mf}=',num2str(round(mean(wmf92*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{['Gross1, R^2 = ',num2str(round(gofGross1.rsquare,3))],...
                    ['Net1, R^2 = ',num2str(round(gofNet1.rsquare,3))],...
                    ['Gross2, R^2 = ',num2str(round(gofGross2.rsquare,3))],...
                    ['Net2, R^2 = ',num2str(round(gofNet2.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';

ax2=axes(fig,'Units','centimeters','Position',[2,5,3,3]);
imshow('withoutBaffle.tiff','Parent',ax2);

saveas(fig,'Figure15.tiff');

save('resultsTRINI.mat','tab92','p_A92','T_A92','wmf92','-append');


%% Finned tubes, fin pitch 9 mm, fin thickness 1 mm (9/1)
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_Finned91_'));


%Initialize table
tab91=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab91.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab91{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:17}=tabloc{:,4:17}+273.15;
    tabloc{:,29:30}=deg2rad(tabloc{:,29:30});
    
    tab91{i,2:30}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab91.p1-tab91.p2+p_amb;
pMain=tab91.p1-tab91.p3+p_amb;
pOut=tab91.p1-tab91.p4+p_amb;

rhoIn=pIn./(R_A.*(tab91.T2));
rhoMain=pMain./(R_A.*(tab91.T3));
rhoOut=pOut./(R_A.*(tab91.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab91.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab91.p6+tab91.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab91.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab91.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

p_A91=pMain-mean([tab91.p6,tab91.p7],2);
TAin=tab91.T3;
T_A91=mean([TAin,tab91.TAout],2);
rho_A=DryAir.rho(p_A91,T_A91);
wmf91=w_mf(d_p,rho_p,p_A91,T_A91);
Vdot_mfMain=wmf91.*A_FBmain;
tab91.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FG91=tab91.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tab91.T_bed11,tab91.T_bed12,tab91.T_bed13,tab91.T_bed14],2);
T_bed2=mean([tab91.T_bed21,tab91.T_bed22,tab91.T_bed23,tab91.T_bed24],2);

tab91.alpha_gross1=tab91.U1.*tab91.I1.*cos(tab91.phi1)./A./(tab91.T_surf1-T_bed1);
tab91.alpha_gross2=tab91.U2.*tab91.I2.*cos(tab91.phi2)./A./(tab91.T_surf2-T_bed2);
tab91.alpha_grossMean=mean([tab91.alpha_gross1,tab91.alpha_gross2],2);
alpha_grossMean91=tab91.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab91.TAout)-DryAir.h(TAin));
tab91.alpha_net1=(tab91.U1.*tab91.I1.*cos(tab91.phi1)-QdotLoss./2)./A./(tab91.T_surf1-T_bed1);
tab91.alpha_net2=(tab91.U2.*tab91.I2.*cos(tab91.phi2)-QdotLoss./2)./A./(tab91.T_surf2-T_bed2);
tab91.alpha_netMean=mean([tab91.alpha_net1,tab91.alpha_net2],2);
alpha_netMean91=tab91.alpha_netMean;


%Calculate accuracy
tab91=accTRINI(tab91,false);


%Create figure
fig=figure(16);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[2,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=tab91.FG;
FGbounds=[2,2.6,3.3,4,5,6];
[gofGross1,gofNet1]=plotGrossNet(ax,FG,FGbounds,tab91.alpha_gross1,tab91.alpha_net1,xfit,[1,2]);
[gofGross2,gofNet2]=plotGrossNet(ax,FG,FGbounds,tab91.alpha_gross2,tab91.alpha_net2,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['Virtual HTCs TRINI Rig, Finned 9/1, w_{mf}=',num2str(round(mean(wmf91*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{['Gross1, R^2 = ',num2str(round(gofGross1.rsquare,3))],...
                    ['Net1, R^2 = ',num2str(round(gofNet1.rsquare,3))],...
                    ['Gross2, R^2 = ',num2str(round(gofGross2.rsquare,3))],...
                    ['Net2, R^2 = ',num2str(round(gofNet2.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';

ax2=axes(fig,'Units','centimeters','Position',[1.8,0.6,3,3]);
imshow('withoutBaffle.tiff','Parent',ax2);

saveas(fig,'Figure16.tiff');

save('resultsTRINI.mat','tab91','p_A91','T_A91','wmf91','-append');


%% Finned tubes, comparison between 9/2 and 9/1
%Create figure
fig=figure(17);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[2,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=[tab92.FG;tab92.FG];
alpha_gross=[tab92.alpha_gross1;tab92.alpha_gross2];
alpha_net=[tab92.alpha_net1;tab92.alpha_net2];
FGbounds=[2,2.5,3,3.7,4.5,6];
[gofGross92,gofNet92]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,[1,2]);


FG=[tab91.FG;tab91.FG];
alpha_gross=[tab91.alpha_gross1;tab91.alpha_gross2];
alpha_net=[tab91.alpha_net1;tab91.alpha_net2];
FGbounds=[2,2.6,3.3,4,5,6];
[gofGross91,gofNet91]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'Virtual HTCs TRINI Rig, Finned 9/2 and 9/1');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(ax,legItems,{['Gross 9/2, R^2 = ',num2str(round(gofGross92.rsquare,3))],...
                    ['Net 9/2, R^2 = ',num2str(round(gofNet92.rsquare,3))],...
                    ['Gross 9/1, R^2 = ',num2str(round(gofGross91.rsquare,3))],...
                    ['Net 9/1, R^2 = ',num2str(round(gofNet91.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure17.tiff');


%% Finned tubes, fin pitch 6 mm, fin thickness 1 mm (6/1)
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_Finned61_'));


%Initialize table
tab61=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab61.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab61{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:17}=tabloc{:,4:17}+273.15;
    tabloc{:,29:30}=deg2rad(tabloc{:,29:30});
    
    tab61{i,2:30}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab61.p1-tab61.p2+p_amb;
pMain=tab61.p1-tab61.p3+p_amb;
pOut=tab61.p1-tab61.p4+p_amb;

rhoIn=pIn./(R_A.*(tab61.T2));
rhoMain=pMain./(R_A.*(tab61.T3));
rhoOut=pOut./(R_A.*(tab61.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab61.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab61.p6+tab61.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab61.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab61.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

p_A=pMain-mean([tab61.p6,tab61.p7],2);
TAin=tab61.T3;
T_A=mean([TAin,tab61.TAout],2);
rho_A=DryAir.rho(p_A,T_A);
wmf=w_mf(d_p,rho_p,p_A,T_A);
Vdot_mfMain=wmf.*A_FBmain;
tab61.FG=mDotAmain./(rho_A.*Vdot_mfMain);


%Calculate gross heat transfer coefficient
T_bed1=mean([tab61.T_bed11,tab61.T_bed12,tab61.T_bed13,tab61.T_bed14],2);
T_bed2=mean([tab61.T_bed21,tab61.T_bed22,tab61.T_bed23,tab61.T_bed24],2);

tab61.alpha_gross1=tab61.U1.*tab61.I1.*cos(tab61.phi1)./A./(tab61.T_surf1-T_bed1);
tab61.alpha_gross2=tab61.U2.*tab61.I2.*cos(tab61.phi2)./A./(tab61.T_surf2-T_bed2);
tab61.alpha_grossMean=mean([tab61.alpha_gross1,tab61.alpha_gross2],2);


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab61.TAout)-DryAir.h(TAin));
tab61.alpha_net1=(tab61.U1.*tab61.I1.*cos(tab61.phi1)-QdotLoss./2)./A./(tab61.T_surf1-T_bed1);
tab61.alpha_net2=(tab61.U2.*tab61.I2.*cos(tab61.phi2)-QdotLoss./2)./A./(tab61.T_surf2-T_bed2);
tab61.alpha_netMean=mean([tab61.alpha_net1,tab61.alpha_net2],2);


%Calculate accuracy
tab61=accTRINI(tab61,false);


%Create figure
fig=figure(18);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[2,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=tab61.FG;
FGbounds=[2,2.7,3.2,3.9,4.5,5.6,6];
[gofGross1,gofNet1]=plotGrossNet(ax,FG,FGbounds,tab61.alpha_gross1,tab61.alpha_net1,xfit,[1,2],'poly1');
[gofGross2,gofNet2]=plotGrossNet(ax,FG,FGbounds,tab61.alpha_gross2,tab61.alpha_net2,xfit,[3,4],'poly1');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['Virtual HTCs TRINI Rig, Finned 6/1, w_{mf}=',num2str(round(mean(wmf*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{['Gross1, R^2 = ',num2str(round(gofGross1.rsquare,3))],...
                    ['Net1, R^2 = ',num2str(round(gofNet1.rsquare,3))],...
                    ['Gross2, R^2 = ',num2str(round(gofGross2.rsquare,3))],...
                    ['Net2, R^2 = ',num2str(round(gofNet2.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';

ax2=axes(fig,'Units','centimeters','Position',[5,5,3,3]);
imshow('withoutBaffle.tiff','Parent',ax2);

saveas(fig,'Figure18.tiff');


%% Finned tubes, comparison between 9/2 and plain tube
%Create figure
fig=figure(19);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6.5];
xfit=linspace(xlim(1),xlim(2),50);
hold on


FG=[tab92.FG;tab92.FG];
alpha_gross=[tab92.alpha_gross1;tab92.alpha_gross2];
alpha_net=[tab92.alpha_net1;tab92.alpha_net2];
FGbounds=[2,2.5,3,3.7,4.5,6];
[gofGross92,gofNet92]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,[1,2]);

FG=[tabPlainWithout.FG;tabPlainWithout.FG];
alpha_gross=[tabPlainWithout.alpha_gross1;tabPlainWithout.alpha_gross2];
alpha_net=[tabPlainWithout.alpha_net1;tabPlainWithout.alpha_net2];
FGbounds=[1.5,2,2.5,3.3,3.6,4,4.5,5.2,5.7,6.5];
[gofGrossWithout,gofNetWithout]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,[3,4]);


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'(Virtual) HTCs TRINI Rig, Plain and Finned 9/2');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'(Virtual) HTC (W/m²K)');
legend(ax,legItems,{['Gross 9/2, R^2 = ',num2str(round(gofGross92.rsquare,3))],...
                    ['Net 9/2, R^2 = ',num2str(round(gofNet92.rsquare,3))],...
                    ['Gross Plain, R^2 = ',num2str(round(gofGrossWithout.rsquare,3))],...
                    ['Net Plain, R^2 = ',num2str(round(gofNetWithout.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure19.tiff');


%% Analyze accuracies
tab=[tabPlainWithout;tabPlainWith;tab92;tab91;tab61];
names={'TRINI_PlainWithout','TRINI_PlainWith','TRINI_Finned92','TRINI_Finned91','TRINI_Finned61'};
longNames={'Plain without Baffle','Plain with Baffle','Finned 9/2','Finned 9/1','Finned 6/1'};


%Initialize table
varnames={'Dataset','FGmin','FGmean','FGmax',...
            'alpha_grossMin','alpha_grossMean','alpha_grossMax',...
            'alpha_netMin','alpha_netMean','alpha_netMax'};
tabAcc=table('Size',[length(names),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tabAcc.Properties.VariableNames=varnames;
tabAcc.Dataset=names';


%Calculate boundaries
tab.FGlower=tab.FGlower./tab.FG;
tab.FGupper=tab.FGupper./tab.FG;

tab.alpha_grossMeanLower=tab.alpha_grossMeanLower./tab.alpha_grossMean;
tab.alpha_grossMeanUpper=tab.alpha_grossMeanUpper./tab.alpha_grossMean;

tab.alpha_netMeanLower=tab.alpha_netMeanLower./tab.alpha_netMean;
tab.alpha_netMeanUpper=tab.alpha_netMeanUpper./tab.alpha_netMean;


%Write table and create figures
for i=1:length(names)
    ind=contains(tab.Dataset,names{i});
    
    tabAcc.FGmin(i)=min([1-tab.FGlower(ind),tab.FGupper(ind)-1],[],'all')*100;
    tabAcc.alpha_grossMin(i)=min([1-tab.alpha_grossMeanLower(ind),tab.alpha_grossMeanUpper(ind)-1],[],'all')*100;
    tabAcc.alpha_netMin(i)=min([1-tab.alpha_netMeanLower(ind),tab.alpha_netMeanUpper(ind)-1],[],'all')*100;
    
    tabAcc.FGmean(i)=mean([1-tab.FGlower(ind),tab.FGupper(ind)-1],'all')*100;
    tabAcc.alpha_grossMean(i)=mean([1-tab.alpha_grossMeanLower(ind),tab.alpha_grossMeanUpper(ind)-1],'all')*100;
    tabAcc.alpha_netMean(i)=mean([1-tab.alpha_netMeanLower(ind),tab.alpha_netMeanUpper(ind)-1],'all')*100;
    
    tabAcc.FGmax(i)=max([1-tab.FGlower(ind),tab.FGupper(ind)-1],[],'all')*100;
    tabAcc.alpha_grossMax(i)=max([1-tab.alpha_grossMeanLower(ind),tab.alpha_grossMeanUpper(ind)-1],[],'all')*100;
    tabAcc.alpha_netMax(i)=max([1-tab.alpha_netMeanLower(ind),tab.alpha_netMeanUpper(ind)-1],[],'all')*100;
    
    
    %Create figure
    fig=figure(300+i);
    clf(fig);
    ax=gca;
    colors=ax.ColorOrder;
    xlim=[1,nnz(ind)];
    hold on
    
    plot([tab.FGlower(ind),tab.FGupper(ind)],'Color',colors(1,:),'LineStyle','-');
    plot([tab.alpha_grossMeanLower(ind),tab.alpha_grossMeanUpper(ind)],'Color',colors(2,:),'LineStyle','--');
    plot([tab.alpha_netMeanLower(ind),tab.alpha_netMeanUpper(ind)],'Color',colors(3,:),'LineStyle',':','LineWidth',1);
    
    
    %Figure formatting
    legItems=repmat(line(),3,1);
    legItems(1)=plot(NaN,'Color',colors(1,:),'LineStyle','-');
    legItems(2)=plot(NaN,'Color',colors(2,:),'LineStyle','--');
    legItems(3)=plot(NaN,'Color',colors(3,:),'LineStyle',':','LineWidth',1);
    hold off
    
    title(ax,['Measurement Uncertainties, TRINI, ',longNames{i}]);
    xlabel(ax,'Dataset Index (-)');
    ylabel(ax,'Relative Uncertainty (-)');
    legend(legItems,{'FG','\alpha_{gross,mean}','\alpha_{net,mean}'},'Location','bestoutside');
    
    fig.Units='centimeters';
    fig.Position=[10,5,17,8.5];
    ax.XLim=xlim;
    ax.YGrid='on';
    saveas(fig,['Figure',num2str(300+i),'.tiff']);  
end


%% Auxiliary functions
function [gofGross,gofNet]=plotGrossNet(ax,FG,FGbounds,alpha_gross,alpha_net,xfit,idx,fitType)
    if nargin<8
        fitType='poly2';
    end
    
    marker={'o','+','x','s'};
    markersz=[sqrt(15),sqrt(15),sqrt(40),sqrt(30)];
    linestyle={'-','--',':','-.'};
    linewidth=[0.5,0.5,1,0.5];
    colors=ax.ColorOrder;
    
    
    FGcats=FG<FGbounds(2:end) & FG>FGbounds(1:end-1);

    
    %FG
    FGmean=arrayfun(@(x) mean(FG(FGcats(:,x))),1:size(FGcats,2));
    FGpos=arrayfun(@(x) max(FG(FGcats(:,x))),1:size(FGcats,2))-FGmean;
    FGneg=FGmean-arrayfun(@(x) min(FG(FGcats(:,x))),1:size(FGcats,2));
    

    %Gross
    alpha=arrayfun(@(x) mean(alpha_gross(FGcats(:,x))),1:size(FGcats,2));
    alphapos=arrayfun(@(x) max(alpha_gross(FGcats(:,x))),1:size(FGcats,2))-alpha;
    alphaneg=alpha-arrayfun(@(x) min(alpha_gross(FGcats(:,x))),1:size(FGcats,2));

    errorbar(ax,FGmean,alpha,alphaneg,alphapos,FGneg,FGpos,...
                'Color',colors(idx(1),:),'CapSize',3,'LineStyle','none','Marker',marker{idx(1)},'MarkerSize',markersz(idx(1)));

    [fitGross,gofGross]=fit(FG,alpha_gross,fitType);
    plot(ax,xfit,fitGross(xfit),'Color',colors(idx(1),:),'LineStyle',linestyle{idx(1)},'LineWidth',linewidth(idx(1)));
    
    
    %Net
    alpha=arrayfun(@(x) mean(alpha_net(FGcats(:,x))),1:size(FGcats,2));
    alphapos=arrayfun(@(x) max(alpha_net(FGcats(:,x))),1:size(FGcats,2))-alpha;
    alphaneg=alpha-arrayfun(@(x) min(alpha_net(FGcats(:,x))),1:size(FGcats,2));

    errorbar(ax,FGmean,alpha,alphaneg,alphapos,FGneg,FGpos,...
                'Color',colors(idx(2),:),'CapSize',3,'LineStyle','none','Marker',marker{idx(2)},'MarkerSize',markersz(idx(2)));

    [fitNet,gofNet]=fit(FG,alpha_net,fitType);
    plot(ax,xfit,fitNet(xfit),'Color',colors(idx(2),:),'LineStyle',linestyle{idx(2)},'LineWidth',linewidth(idx(2)));
end





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
%https://doi.org/10.5281/zenodo.5474020
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5500329
%
%
%This script analyzes the TRINI test rig data and creates all published
%figures
%Requires all TRINI data files in the same folder and all auxiliary
%functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
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

rho_Anorm=1.293;
eta_A=18.107811e-6;
R_A=287.0533;

p_amb=101325;

alpha=10e-12;
beta=83e-12;
s=20e-3;


dirCont=dir();


%% Plain tubes, with baffle
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_PlainWith_'));


%Initialize table
varnames={'Dataset','P_el1','P_el2','T2','T3',...
            'T4','TAout','T_bed11','T_bed12','T_bed21',...
            'T_bed22','T_surf1','T_surf2','mDotA','p1',...
            'p2','p3','p4','p5','p6',...
            'p7','p8',...
            'FG','alpha_gross1','alpha_gross2','alpha_grossMean','alpha_net1',...
            'alpha_net2','alpha_netMean'};
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:13}=tabloc{:,4:13}+273.15;
    
    tab{i,2:22}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab.p1+tab.p2+p_amb;
pMain=tab.p1+tab.p3+p_amb;
pOut=tab.p1+tab.p4+p_amb;

rhoIn=pIn./(R_A.*(tab.T2));
rhoMain=pMain./(R_A.*(tab.T3));
rhoOut=pOut./(R_A.*(tab.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab.p6+tab.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

TAin=tab.T3;
T_A=mean([TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot_mfMain=w_mf(d_p,rho_p,T_A).*A_FBmain;
tab.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FGwith=tab.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tab.T_bed11,tab.T_bed12],2);
T_bed2=mean([tab.T_bed21,tab.T_bed22],2);

tab.alpha_gross1=tab.P_el1./A./(tab.T_surf1-T_bed1);
tab.alpha_gross2=tab.P_el2./A./(tab.T_surf2-T_bed2);
tab.alpha_grossMean=mean([tab.alpha_gross1,tab.alpha_gross2],2);
alpha_grossMeanWith=tab.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab.TAout)-DryAir.h(TAin));
tab.alpha_net1=(tab.P_el1-QdotLoss./2)./A./(tab.T_surf1-T_bed1);
tab.alpha_net2=(tab.P_el2-QdotLoss./2)./A./(tab.T_surf2-T_bed2);
tab.alpha_netMean=mean([tab.alpha_net1,tab.alpha_net2],2);
alpha_netMeanWith=tab.alpha_netMean;


%Create figure
fig=figure(14);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross1;
scatter(ax,x,y,10,colors(1,:),'o');

fitGross1=fit(x,y,'poly2');
plot(xfit,fitGross1(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net1;
scatter(ax,x,y,10,colors(2,:),'+');

fitNet1=fit(x,y,'poly2');
plot(xfit,fitNet1(xfit),'Color',colors(2,:),'LineStyle','--');


y=tab.alpha_gross2;
scatter(ax,x,y,20,colors(3,:),'x');

fitGross2=fit(x,y,'poly2');
plot(xfit,fitGross2(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net2;
scatter(ax,x,y,20,colors(4,:),'s');

fitNet2=fit(x,y,'poly2');
plot(xfit,fitNet2(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'HTCs TRINI Rig, Plain with Baffle');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{'Gross1','Net1','Gross2','Net2'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure14.tiff');


%% Plain tubes, without baffle
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_PlainWithout_'));


%Initialize table
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:13}=tabloc{:,4:13}+273.15;
    
    tab{i,2:22}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab.p1+tab.p2+p_amb;
pMain=tab.p1+tab.p3+p_amb;
pOut=tab.p1+tab.p4+p_amb;

rhoIn=pIn./(R_A.*(tab.T2));
rhoMain=pMain./(R_A.*(tab.T3));
rhoOut=pOut./(R_A.*(tab.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab.p6+tab.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

TAin=tab.T3;
T_A=mean([TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot_mfMain=w_mf(d_p,rho_p,T_A).*A_FBmain;
tab.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FGwithout=tab.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tab.T_bed11,tab.T_bed12],2);
T_bed2=mean([tab.T_bed21,tab.T_bed22],2);

tab.alpha_gross1=tab.P_el1./A./(tab.T_surf1-T_bed1);
tab.alpha_gross2=tab.P_el2./A./(tab.T_surf2-T_bed2);
tab.alpha_grossMean=mean([tab.alpha_gross1,tab.alpha_gross2],2);
alpha_grossMeanWithout=tab.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab.TAout)-DryAir.h(TAin));
tab.alpha_net1=(tab.P_el1-QdotLoss./2)./A./(tab.T_surf1-T_bed1);
tab.alpha_net2=(tab.P_el2-QdotLoss./2)./A./(tab.T_surf2-T_bed2);
tab.alpha_netMean=mean([tab.alpha_net1,tab.alpha_net2],2);
alpha_netMeanWithout=tab.alpha_netMean;


%Create figure
fig=figure(15);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross1;
scatter(ax,x,y,10,colors(1,:),'o');

fitGross1=fit(x,y,'poly2');
plot(xfit,fitGross1(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net1;
scatter(ax,x,y,10,colors(2,:),'+');

fitNet1=fit(x,y,'poly2');
plot(xfit,fitNet1(xfit),'Color',colors(2,:),'LineStyle','--');


y=tab.alpha_gross2;
scatter(ax,x,y,20,colors(3,:),'x');

fitGross2=fit(x,y,'poly2');
plot(xfit,fitGross2(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net2;
scatter(ax,x,y,20,colors(4,:),'s');

fitNet2=fit(x,y,'poly2');
plot(xfit,fitNet2(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'HTCs TRINI Rig, Plain without Baffle');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{'Gross1','Net1','Gross2','Net2'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure15.tiff');


%% Plain tubes, comparison with / without baffle
%Create figure
fig=figure(16);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
xfit=linspace(xlim(1),xlim(2),50);
hold on


x=FGwith;
y=alpha_grossMeanWith;
fitGrossWith=fit(x,y,'poly2');
plot(xfit,fitGrossWith(xfit),'Color',colors(1,:),'LineStyle','-');


y=alpha_netMeanWith;
fitNetWith=fit(x,y,'poly2');
plot(xfit,fitNetWith(xfit),'Color',colors(2,:),'LineStyle','--');


x=FGwithout;
y=alpha_grossMeanWithout;
fitGrossWithout=fit(x,y,'poly2');
plot(xfit,fitGrossWithout(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=alpha_netMeanWithout;
fitNetWithout=fit(x,y,'poly2');
plot(xfit,fitNetWithout(xfit),'Color',colors(4,:),'LineStyle','-.');
hold off


%Figure formatting
title(ax,'HTCs TRINI Rig, Plain with and without Baffle');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,{'Gross With','Net With','Gross Without','Net Without'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure16.tiff');


%% Finned tubes, fin pitch 9 mm, fin thickness 2 mm (9/2)
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_Finned92_'));


%Initialize table
varnames={'Dataset','P_el1','P_el2','T2','T3',...
            'T4','TAout','T_bed11','T_bed12','T_bed13',...
            'T_bed14','T_bed21','T_bed22','T_bed23','T_bed24',...
            'T_surf1','T_surf2','mDotA','p1',...
            'p2','p3','p4','p5','p6',...
            'p7','p8',...
            'FG','alpha_gross1','alpha_gross2','alpha_grossMean','alpha_net1',...
            'alpha_net2','alpha_netMean'};
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:17}=tabloc{:,4:17}+273.15;
    
    tab{i,2:26}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab.p1+tab.p2+p_amb;
pMain=tab.p1+tab.p3+p_amb;
pOut=tab.p1+tab.p4+p_amb;

rhoIn=pIn./(R_A.*(tab.T2));
rhoMain=pMain./(R_A.*(tab.T3));
rhoOut=pOut./(R_A.*(tab.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab.p6+tab.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

TAin=tab.T3;
T_A=mean([TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot_mfMain=w_mf(d_p,rho_p,T_A).*A_FBmain;
tab.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FG92=tab.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tab.T_bed11,tab.T_bed12,tab.T_bed13,tab.T_bed14],2);
T_bed2=mean([tab.T_bed21,tab.T_bed22,tab.T_bed23,tab.T_bed24],2);

tab.alpha_gross1=tab.P_el1./A./(tab.T_surf1-T_bed1);
tab.alpha_gross2=tab.P_el2./A./(tab.T_surf2-T_bed2);
tab.alpha_grossMean=mean([tab.alpha_gross1,tab.alpha_gross2],2);
alpha_grossMean92=tab.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab.TAout)-DryAir.h(TAin));
tab.alpha_net1=(tab.P_el1-QdotLoss./2)./A./(tab.T_surf1-T_bed1);
tab.alpha_net2=(tab.P_el2-QdotLoss./2)./A./(tab.T_surf2-T_bed2);
tab.alpha_netMean=mean([tab.alpha_net1,tab.alpha_net2],2);
alpha_netMean92=tab.alpha_netMean;


%Create figure
fig=figure(17);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross1;
scatter(ax,x,y,10,colors(1,:),'o');

fitGross1=fit(x,y,'poly2');
plot(xfit,fitGross1(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net1;
scatter(ax,x,y,10,colors(2,:),'+');

fitNet1=fit(x,y,'poly2');
plot(xfit,fitNet1(xfit),'Color',colors(2,:),'LineStyle','--');


y=tab.alpha_gross2;
scatter(ax,x,y,20,colors(3,:),'x');

fitGross2=fit(x,y,'poly2');
plot(xfit,fitGross2(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net2;
scatter(ax,x,y,20,colors(4,:),'s');

fitNet2=fit(x,y,'poly2');
plot(xfit,fitNet2(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'Virtual HTCs TRINI Rig, Finned 9/2');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{'Gross1','Net1','Gross2','Net2'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure17.tiff');


%% Finned tubes, fin pitch 9 mm, fin thickness 1 mm (9/1)
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_Finned91_'));


%Initialize table
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:17}=tabloc{:,4:17}+273.15;
    
    tab{i,2:26}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab.p1+tab.p2+p_amb;
pMain=tab.p1+tab.p3+p_amb;
pOut=tab.p1+tab.p4+p_amb;

rhoIn=pIn./(R_A.*(tab.T2));
rhoMain=pMain./(R_A.*(tab.T3));
rhoOut=pOut./(R_A.*(tab.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab.p6+tab.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

TAin=tab.T3;
T_A=mean([TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot_mfMain=w_mf(d_p,rho_p,T_A).*A_FBmain;
tab.FG=mDotAmain./(rho_A.*Vdot_mfMain);
FG91=tab.FG;


%Calculate gross heat transfer coefficient
T_bed1=mean([tab.T_bed11,tab.T_bed12,tab.T_bed13,tab.T_bed14],2);
T_bed2=mean([tab.T_bed21,tab.T_bed22,tab.T_bed23,tab.T_bed24],2);

tab.alpha_gross1=tab.P_el1./A./(tab.T_surf1-T_bed1);
tab.alpha_gross2=tab.P_el2./A./(tab.T_surf2-T_bed2);
tab.alpha_grossMean=mean([tab.alpha_gross1,tab.alpha_gross2],2);
alpha_grossMean91=tab.alpha_grossMean;


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab.TAout)-DryAir.h(TAin));
tab.alpha_net1=(tab.P_el1-QdotLoss./2)./A./(tab.T_surf1-T_bed1);
tab.alpha_net2=(tab.P_el2-QdotLoss./2)./A./(tab.T_surf2-T_bed2);
tab.alpha_netMean=mean([tab.alpha_net1,tab.alpha_net2],2);
alpha_netMean91=tab.alpha_netMean;


%Create figure
fig=figure(18);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross1;
scatter(ax,x,y,10,colors(1,:),'o');

fitGross1=fit(x,y,'poly2');
plot(xfit,fitGross1(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net1;
scatter(ax,x,y,10,colors(2,:),'+');

fitNet1=fit(x,y,'poly2');
plot(xfit,fitNet1(xfit),'Color',colors(2,:),'LineStyle','--');


y=tab.alpha_gross2;
scatter(ax,x,y,20,colors(3,:),'x');

fitGross2=fit(x,y,'poly2');
plot(xfit,fitGross2(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net2;
scatter(ax,x,y,20,colors(4,:),'s');

fitNet2=fit(x,y,'poly2');
plot(xfit,fitNet2(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'Virtual HTCs TRINI Rig, Finned 9/1');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{'Gross1','Net1','Gross2','Net2'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure18.tiff');


%% Finned tubes, comparison between 9/2 and 9/1
%Create figure
fig=figure(19);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
xfit=linspace(xlim(1),xlim(2),50);
hold on


x=FG92;
y=alpha_grossMean92;
fitGrossWith=fit(x,y,'poly2');
plot(xfit,fitGrossWith(xfit),'Color',colors(1,:),'LineStyle','-');


y=alpha_netMean92;
fitNetWith=fit(x,y,'poly2');
plot(xfit,fitNetWith(xfit),'Color',colors(2,:),'LineStyle','--');


x=FG91;
y=alpha_grossMean91;
fitGrossWithout=fit(x,y,'poly2');
plot(xfit,fitGrossWithout(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=alpha_netMean91;
fitNetWithout=fit(x,y,'poly2');
plot(xfit,fitNetWithout(xfit),'Color',colors(4,:),'LineStyle','-.');
hold off


%Figure formatting
title(ax,'Virtual HTCs TRINI Rig, Finned 9/2 and 9/1');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(ax,{'Gross 9/2','Net 9/2','Gross 9/1','Net 9/1'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure19.tiff');


%% Finned tubes, fin pitch 6 mm, fin thickness 1 mm (6/1)
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'TRINI_Finned61_'));


%Initialize table
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,4:17}=tabloc{:,4:17}+273.15;
    
    tab{i,2:26}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
pIn=tab.p1+tab.p2+p_amb;
pMain=tab.p1+tab.p3+p_amb;
pOut=tab.p1+tab.p4+p_amb;

rhoIn=pIn./(R_A.*(tab.T2));
rhoMain=pMain./(R_A.*(tab.T3));
rhoOut=pOut./(R_A.*(tab.T4));

lambdaIn=eta_A.*beta.*A_FBin./(rhoIn.*alpha);
lambdaMain=eta_A.*beta.*A_FBmain./(rhoMain.*alpha);
lambdaOut=eta_A.*beta.*A_FBout./(rhoOut.*alpha);

VdotIn_est=-lambdaIn./2+sqrt((lambdaIn./2).^2+beta.*A_FBin.^2.*tab.p5./(rhoIn.*s));
VdotMain_est=-lambdaMain./2+sqrt((lambdaMain./2).^2+beta.*A_FBmain.^2.*(tab.p6+tab.p7)./(2*rhoMain.*s));
VdotOut_est=-lambdaOut./2+sqrt((lambdaOut./2).^2+beta.*A_FBout.^2.*tab.p8./(rhoOut.*s));

mDotAin_est=VdotIn_est.*rhoIn;
mDotAmain_est=VdotMain_est.*rhoMain;
mDotAout_est=VdotOut_est.*rhoOut;

mDotA_est=mDotAin_est+mDotAmain_est+mDotAout_est;

% mDotAin=mDotAin_est./mDotA_est.*tab.mDotA;
mDotAmain=mDotAmain_est./mDotA_est.*tab.mDotA;
% mDotAout=mDotAout_est./mDotA_est.*tab.mDotA;

TAin=tab.T3;
T_A=mean([TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot_mfMain=w_mf(d_p,rho_p,T_A).*A_FBmain;
tab.FG=mDotAmain./(rho_A.*Vdot_mfMain);


%Calculate gross heat transfer coefficient
T_bed1=mean([tab.T_bed11,tab.T_bed12,tab.T_bed13,tab.T_bed14],2);
T_bed2=mean([tab.T_bed21,tab.T_bed22,tab.T_bed23,tab.T_bed24],2);

tab.alpha_gross1=tab.P_el1./A./(tab.T_surf1-T_bed1);
tab.alpha_gross2=tab.P_el2./A./(tab.T_surf2-T_bed2);
tab.alpha_grossMean=mean([tab.alpha_gross1,tab.alpha_gross2],2);


%Calculate net heat transfer coefficient
QdotLoss=mDotAmain.*(DryAir.h(tab.TAout)-DryAir.h(TAin));
tab.alpha_net1=(tab.P_el1-QdotLoss./2)./A./(tab.T_surf1-T_bed1);
tab.alpha_net2=(tab.P_el2-QdotLoss./2)./A./(tab.T_surf2-T_bed2);
tab.alpha_netMean=mean([tab.alpha_net1,tab.alpha_net2],2);


%Create figure
fig=figure(20);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross1;
scatter(ax,x,y,10,colors(1,:),'o');

fitGross1=fit(x,y,'poly1');
plot(xfit,fitGross1(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net1;
scatter(ax,x,y,10,colors(2,:),'+');

fitNet1=fit(x,y,'poly1');
plot(xfit,fitNet1(xfit),'Color',colors(2,:),'LineStyle','--');


y=tab.alpha_gross2;
scatter(ax,x,y,20,colors(3,:),'x');

fitGross2=fit(x,y,'poly1');
plot(xfit,fitGross2(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net2;
scatter(ax,x,y,20,colors(4,:),'s');

fitNet2=fit(x,y,'poly1');
plot(xfit,fitNet2(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'Virtual HTCs TRINI Rig, Finned 6/1');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{'Gross1','Net1','Gross2','Net2'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure20.tiff');


%% Finned tubes, comparison between 9/2 and 9/1
%Create figure
fig=figure(21);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,7];
xfit=linspace(xlim(1),xlim(2),50);
hold on


x=FG92;
y=alpha_grossMean92;
fitGrossWith=fit(x,y,'poly2');
plot(xfit,fitGrossWith(xfit),'Color',colors(1,:),'LineStyle','-');


y=alpha_netMean92;
fitNetWith=fit(x,y,'poly2');
plot(xfit,fitNetWith(xfit),'Color',colors(2,:),'LineStyle','--');


x=FGwithout;
y=alpha_grossMeanWithout;
fitGrossWithout=fit(x,y,'poly2');
plot(xfit,fitGrossWithout(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=alpha_netMeanWithout;
fitNetWithout=fit(x,y,'poly2');
plot(xfit,fitNetWithout(xfit),'Color',colors(4,:),'LineStyle','-.');
hold off


%Figure formatting
title(ax,'(Virtual) HTCs TRINI Rig, Plain and Finned 9/2');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'(Virtual) HTC (W/m²K)');
legend(ax,{'Gross 9/2','Net 9/2','Gross Plain','Net Plain'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure21.tiff');





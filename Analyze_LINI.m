%% Analyze LINI test rig results
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
%This script analyzes the LINI test rig data and creates all published
%figures
%Requires all LINI data files in the same folder and all auxiliary
%functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
%   - w_mf.m


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

d_p=146e-6;
rho_p=2650;
rho_Anorm=1.293;


dirCont=dir();


%% Plain tubes
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'LINI_Plain_'));


%Initialize table
varnames={'Dataset','P_el','TAin','TAout','T_bed',...
            'T_surf','VdotRota','pRota','Vdot',...
            'FG','alpha_gross','alpha_net'};
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,3:6}=tabloc{:,3:6}+273.15;
    tabloc{:,[7,9]}=tabloc{:,[7,9]}./60^2;
    
    tab{i,2:9}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
T_A=mean([tab.TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot=tab.Vdot.*rho_Anorm./rho_A;
Vdot_mf=w_mf(d_p,rho_p,T_A).*A_FB;
tab.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
tab.alpha_gross=tab.P_el./A./(tab.T_surf-tab.T_bed);


%Calculate net heat transfer coefficient
mDotA_loss=Vdot.*rho_A.*x_heated;
QdotLoss=mDotA_loss.*(DryAir.h(tab.TAout)-DryAir.h(tab.TAin));
tab.alpha_net=(tab.P_el-QdotLoss)./A./(tab.T_surf-tab.T_bed);


%Create figure
fig=figure(11);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross;
scatter(ax,x,y,10,colors(1,:),'o');

fitPlainGross=fit(x,y,'poly2');
plot(xfit,fitPlainGross(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net;
scatter(ax,x,y,10,colors(2,:),'+');

fitPlainNet=fit(x,y,'poly2');
plot(xfit,fitPlainNet(xfit),'Color',colors(2,:),'LineStyle','--');


%Figure formatting
legItems=repmat(line(),2,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
hold off

title(ax,'HTCs LINI Rig, Plain Tube');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{'Gross','Net'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure11.tiff');


%% Finned tubes
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'LINI_Finned_'));


%Initialize table
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,3:6}=tabloc{:,3:6}+273.15;
    tabloc{:,[7,9]}=tabloc{:,[7,9]}./60^2;
    
    tab{i,2:9}=mean(tabloc{:,2:end});
end


%Calculate degree of fluidization
T_A=mean([tab.TAin,tab.TAout],2);
rho_A=DryAir.rho(T_A);
Vdot=tab.Vdot.*rho_Anorm./rho_A;
Vdot_mf=w_mf(d_p,rho_p,T_A).*A_FB;
tab.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
tab.alpha_gross=tab.P_el./A./(tab.T_surf-tab.T_bed);


%Calculate net heat transfer coefficient
mDotA_loss=Vdot.*rho_A.*x_heated;
QdotLoss=mDotA_loss.*(DryAir.h(tab.TAout)-DryAir.h(tab.TAin));
tab.alpha_net=(tab.P_el-QdotLoss)./A./(tab.T_surf-tab.T_bed);


%Create figure
fig=figure(12);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6];
x=tab.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tab.alpha_gross;
scatter(ax,x,y,10,colors(1,:),'o');

fitFinnedGross=fit(x,y,'poly2');
plot(xfit,fitFinnedGross(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net;
scatter(ax,x,y,10,colors(2,:),'+');

fitFinnedNet=fit(x,y,'poly2');
plot(xfit,fitFinnedNet(xfit),'Color',colors(2,:),'LineStyle','--');


%Figure formatting
legItems=repmat(line(),2,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
hold off

title(ax,'Virtual HTCs LINI Rig, Finned Tube');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{'Gross','Net'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure12.tiff');


%% Comparison between plain and finned tubes
%Create figure
fig=figure(13);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


%Plot previous fits
plot(xfit,fitPlainGross(xfit),'Color',colors(1,:),'LineStyle','-');
plot(xfit,fitPlainNet(xfit),'Color',colors(2,:),'LineStyle','--');
plot(xfit,fitFinnedGross(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);
plot(xfit,fitFinnedNet(xfit),'Color',colors(4,:),'LineStyle','-.');
hold off


%Figure formatting
title(ax,'(Virtual) HTCs LINI Rig');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'(Virtual) HTC (W/m²K)');
legend(ax,{'Plain Gross','Plain Net','Finned Gross','Finned Net'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure13.tiff');





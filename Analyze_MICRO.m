%% Analyze MICRO test rig results
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
%This script analyzes the MICRO test rig data and creates all published
%figures
%Requires all MICRO data files in the same folder and all auxiliary
%functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - MICRO_MeanValues.csv
%   - @DryAir
%   - w_mf.m


%% Constants
l=0.2;
w=0.2;

l_heated=0.1;
d=25e-3;
A=d*pi*l_heated;

d_p1=87.141e-6;
d_p2=210e-6;

rho_p=2650;
rho_Anorm=1.293;

dirCont=dir();

%Read mean fluidization values (Rotameter)
tabmean=readtable('MICRO_MeanValues.csv','VariableNamingRule','preserve');


%% 87 µm
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'MICRO_87µm_'));


%Initialize table
varnames={'Dataset','T_Lein','T_Laus','T_BED1','T_PROBE1',...
            'T_PROBE2','P_EL','VdotAir','VdotAir15','VdotAir20',...
            'FG','alpha_gross','alpha_net'};
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,2:6}=tabloc{:,2:6}+273.15;
    tab{i,2:7}=mean(tabloc{:,2:end});
    
    tab{i,8:10}=tabmean{matches(tabmean.Dataset,tab{i,1}),2:end}./60^2;
end


%Calculate degree of fluidization
rho_A=DryAir.rho(tab.T_BED1);
Vdot=mean([tab.VdotAir15,tab.VdotAir20],2).*rho_Anorm./rho_A;
Vdot_mf=w_mf(d_p1,rho_p,tab.T_BED1).*l.*w;
tab.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
T_surf=mean([tab.T_PROBE1,tab.T_PROBE2],2);
tab.alpha_gross=tab.P_EL./A./(T_surf-tab.T_BED1);


%Calculate net heat transfer coefficient
mDotA=Vdot.*rho_A;
Qdot_loss=mDotA.*(DryAir.h(tab.T_Laus)-DryAir.h(tab.T_Lein));
tab.alpha_net=(tab.P_EL-Qdot_loss)./A./(T_surf-tab.T_BED1);


%Create figure
fig=figure(8);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[3.5,13];
xfit=linspace(xlim(1),xlim(2),50);
hold on


%Regular Pitch
pitch=contains(tab.Dataset,'87µm_RegularPitch');
x=tab.FG(pitch);


y=tab.alpha_gross(pitch);
scatter(ax,x,y,10,colors(1,:),'o');

yfit=fit(x,y,'poly2');
plot(xfit,yfit(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net(pitch);
scatter(ax,x,y,10,colors(2,:),'+');

yfit=fit(x,y,'poly2');
plot(xfit,yfit(xfit),'Color',colors(2,:),'LineStyle','--');


%sandTES Pitch
pitch=contains(tab.Dataset,'87µm_sandTESPitch');
x=tab.FG(pitch);


y=tab.alpha_gross(pitch);
scatter(ax,x,y,20,colors(3,:),'x');

fitST87gross=fit(x,y,'poly2');
plot(xfit,fitST87gross(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net(pitch);
scatter(ax,x,y,20,colors(4,:),'s');

fitST87net=fit(x,y,'poly2');
plot(xfit,fitST87net(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'HTCs MICRO Rig, d_p=87 µm');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{'Regular Spacing Gross','Regular Spacing Net','sandTES Spacing Gross','sandTES Spacing Net'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure8.tiff');


%% 210 µm
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'MICRO_210µm_'));


%Initialize table
tab=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,2:6}=tabloc{:,2:6}+273.15;
    tab{i,2:7}=mean(tabloc{:,2:end});
    
    tab{i,8:10}=tabmean{matches(tabmean.Dataset,tab{i,1}),2:end}./60^2;
end


%Calculate degree of fluidization
rho_A=DryAir.rho(tab.T_BED1);
Vdot=mean([tab.VdotAir15,tab.VdotAir20],2).*rho_Anorm./rho_A;
Vdot_mf=w_mf(d_p2,rho_p,tab.T_BED1).*l.*w;
tab.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
T_surf=mean([tab.T_PROBE1,tab.T_PROBE2],2);
tab.alpha_gross=tab.P_EL./A./(T_surf-tab.T_BED1);


%Calculate net heat transfer coefficient
mDotA=Vdot.*rho_A;
Qdot_loss=mDotA.*(DryAir.h(tab.T_Laus)-DryAir.h(tab.T_Lein));
tab.alpha_net=(tab.P_EL-Qdot_loss)./A./(T_surf-tab.T_BED1);


%Create figure
fig=figure(9);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


%Regular Pitch
pitch=contains(tab.Dataset,'210µm_RegularPitch');
x=tab.FG(pitch);


y=tab.alpha_gross(pitch);
scatter(ax,x,y,10,colors(1,:),'o');

yfit=fit(x,y,'poly2');
plot(xfit,yfit(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab.alpha_net(pitch);
scatter(ax,x,y,10,colors(2,:),'+');

yfit=fit(x,y,'poly2');
plot(xfit,yfit(xfit),'Color',colors(2,:),'LineStyle','--');


%sandTES Pitch
pitch=contains(tab.Dataset,'210µm_sandTESPitch');
x=tab.FG(pitch);


y=tab.alpha_gross(pitch);
scatter(ax,x,y,20,colors(3,:),'x');

fitST210gross=fit(x,y,'poly2');
plot(xfit,fitST210gross(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab.alpha_net(pitch);
scatter(ax,x,y,20,colors(4,:),'s');

fitST210net=fit(x,y,'poly2');
plot(xfit,fitST210net(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,'HTCs MICRO Rig, d_p=210 µm');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{'Regular Spacing Gross','Regular Spacing Net','sandTES Spacing Gross','sandTES Spacing Net'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YLim=[100,350];
ax.YGrid='on';
saveas(fig,'Figure9.tiff');


%% Comparison between 87 µm and 210 µm
%Create figure
fig=figure(10);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[3,7];
xfit=linspace(xlim(1),xlim(2),50);
hold on


%Plot previous fits
plot(xfit,fitST87gross(xfit),'Color',colors(1,:),'LineStyle','-');
plot(xfit,fitST87net(xfit),'Color',colors(2,:),'LineStyle','--');
plot(xfit,fitST210gross(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);
plot(xfit,fitST210net(xfit),'Color',colors(4,:),'LineStyle','-.');
hold off


%Figure formatting
title(ax,'HTCs MICRO Rig, sandTES Spacing');
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,{'87 µm Gross','87 µm Net','210 µm Gross','210 µm Net'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure10.tiff');





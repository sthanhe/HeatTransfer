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
%https://doi.org/10.5281/zenodo.5802409
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5802407
%
%
%
%This script analyzes the MICRO test rig data and creates all published
%figures.
%
%Requires all MICRO data files in the same folder and all auxiliary
%functions on the MATLAB path.
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
%   - accMICRO.m
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
p_amb=1013.25e2;

dirCont=dir();

%Read mean fluidization values (Rotameter)
tabmean=readtable('MICRO_MeanValues.csv','VariableNamingRule','preserve');


%% 87 µm
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'MICRO_87mum_'));


%Initialize table
varnames={'Dataset','T_Lein','T_Laus','T_BED1','T_PROBE1',...
            'T_PROBE2','U','I','P_FEED','VdotAir','VdotAir15',...
            'VdotAir20','FG','FGlower','FGupper',...
                'alpha_gross','alpha_grossLower','alpha_grossUpper',...
                'alpha_net','alpha_netLower','alpha_netUpper'};
tab87=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab87.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab87{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,2:6}=tabloc{:,2:6}+273.15;
    tabloc.P_FEED=tabloc.P_FEED.*10^2;
    tab87{i,2:9}=mean(tabloc{:,2:end});
    
    tab87{i,10:12}=tabmean{matches(tabmean.Dataset,tab87{i,1}),2:end}./60^2;
end


%Calculate degree of fluidization
p_A=p_amb+4335;
rho_A=DryAir.rho(p_A,tab87.T_BED1);
Vdot=mean([tab87.VdotAir15,tab87.VdotAir20],2).*rho_Anorm./rho_A;
wmf87=w_mf(d_p1,rho_p,p_A,tab87.T_BED1);
Vdot_mf=wmf87.*l.*w;
tab87.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
T_surf=mean([tab87.T_PROBE1,tab87.T_PROBE2],2);
tab87.alpha_gross=tab87.U.*tab87.I./A./(T_surf-tab87.T_BED1);


%Calculate net heat transfer coefficient
mDotA=Vdot.*rho_A;
Qdot_loss=mDotA.*(DryAir.h(tab87.T_Laus)-DryAir.h(tab87.T_Lein));
tab87.alpha_net=(tab87.U.*tab87.I-Qdot_loss)./A./(T_surf-tab87.T_BED1);


%Calculate accuracy
tab87=accMICRO(tab87,d_p1);


%Create figure
fig=figure(9);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[3.5,13];
xfit=linspace(xlim(1),xlim(2),50);
hold on


%Regular Pitch
pitch=contains(tab87.Dataset,'87mum_RegularPitch');
x=tab87.FG(pitch);


y=tab87.alpha_gross(pitch);
scatter(ax,x,y,10,colors(1,:),'o');

[fitReg87gross,gofReg87gross]=fit(x,y,'poly2');
plot(xfit,fitReg87gross(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab87.alpha_net(pitch);
scatter(ax,x,y,10,colors(2,:),'+');

[fitReg87net,gofReg87net]=fit(x,y,'poly2');
plot(xfit,fitReg87net(xfit),'Color',colors(2,:),'LineStyle','--');


%sandTES Pitch
pitch=contains(tab87.Dataset,'87mum_sandTESPitch');
x=tab87.FG(pitch);


y=tab87.alpha_gross(pitch);
scatter(ax,x,y,20,colors(3,:),'x');

[fitST87gross,gofST87gross]=fit(x,y,'poly2');
plot(xfit,fitST87gross(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab87.alpha_net(pitch);
scatter(ax,x,y,20,colors(4,:),'s');

[fitST87net,gofST87net]=fit(x,y,'poly2');
plot(xfit,fitST87net(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['HTCs MICRO Rig, d_p=87 µm, w_{mf}=',num2str(round(mean(wmf87*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{['Regular Spacing Gross',newline,'R^2 = ',num2str(round(gofReg87gross.rsquare,3))],...
                    ['Regular Spacing Net',newline,'R^2 = ',num2str(round(gofReg87net.rsquare,3))],...
                    ['sandTES Spacing Gross',newline,'R^2 = ',num2str(round(gofST87gross.rsquare,3))],...
                    ['sandTES Spacing Net',newline,'R^2 = ',num2str(round(gofST87net.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure9.tiff');

save('resultsMICRO.mat','tab87','fitReg87gross','fitST87gross','wmf87');


%% 210 µm
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'MICRO_210mum_'));


%Initialize table
varnames={'Dataset','T_Lein','T_Laus','T_BED1','T_PROBE1',...
            'T_PROBE2','U','I','VdotAir','VdotAir15',...
            'VdotAir20','FG','FGlower','FGupper',...
                'alpha_gross','alpha_grossLower','alpha_grossUpper',...
                'alpha_net','alpha_netLower','alpha_netUpper'};
tab210=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tab210.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tab210{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,2:6}=tabloc{:,2:6}+273.15;
    tab210{i,2:8}=mean(tabloc{:,2:end});
    
    tab210{i,9:11}=tabmean{matches(tabmean.Dataset,tab210{i,1}),2:end}./60^2;
end


%Calculate degree of fluidization
p_A=p_amb+4335;
rho_A=DryAir.rho(p_A,tab210.T_BED1);
Vdot=mean([tab210.VdotAir15,tab210.VdotAir20],2).*rho_Anorm./rho_A;
wmf210=w_mf(d_p2,rho_p,p_A,tab210.T_BED1);
Vdot_mf=wmf210.*l.*w;
tab210.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
T_surf=mean([tab210.T_PROBE1,tab210.T_PROBE2],2);
tab210.alpha_gross=tab210.U.*tab210.I./A./(T_surf-tab210.T_BED1);


%Calculate net heat transfer coefficient
mDotA=Vdot.*rho_A;
Qdot_loss=mDotA.*(DryAir.h(tab210.T_Laus)-DryAir.h(tab210.T_Lein));
tab210.alpha_net=(tab210.U.*tab210.I-Qdot_loss)./A./(T_surf-tab210.T_BED1);


%Calculate accuracy
tab210=accMICRO(tab210,d_p2);


%Create figure
fig=figure(10);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1,6];
xfit=linspace(xlim(1),xlim(2),50);
hold on


%Regular Pitch
pitch=contains(tab210.Dataset,'210mum_RegularPitch');
x=tab210.FG(pitch);


y=tab210.alpha_gross(pitch);
scatter(ax,x,y,10,colors(1,:),'o');

[fitReg210gross,gofReg210gross]=fit(x,y,'poly2');
plot(xfit,fitReg210gross(xfit),'Color',colors(1,:),'LineStyle','-');


y=tab210.alpha_net(pitch);
scatter(ax,x,y,10,colors(2,:),'+');

[fitReg210net,gofReg210net]=fit(x,y,'poly2');
plot(xfit,fitReg210net(xfit),'Color',colors(2,:),'LineStyle','--');


%sandTES Pitch
pitch=contains(tab210.Dataset,'210mum_sandTESPitch');
x=tab210.FG(pitch);


y=tab210.alpha_gross(pitch);
scatter(ax,x,y,20,colors(3,:),'x');

[fitST210gross,gofST210gross]=fit(x,y,'poly2');
plot(xfit,fitST210gross(xfit),'Color',colors(3,:),'LineStyle',':','LineWidth',1);


y=tab210.alpha_net(pitch);
scatter(ax,x,y,20,colors(4,:),'s');

[fitST210net,gofST210net]=fit(x,y,'poly2');
plot(xfit,fitST210net(xfit),'Color',colors(4,:),'LineStyle','-.');


%Figure formatting
legItems=repmat(line(),4,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
legItems(3)=plot(NaN,NaN,'Color',colors(3,:),'LineStyle',':','Marker','x','LineWidth',1);
legItems(4)=plot(NaN,NaN,'Color',colors(4,:),'LineStyle','-.','Marker','s');
hold off

title(ax,['HTCs MICRO Rig, d_p=210 µm, w_{mf}=',num2str(round(mean(wmf210*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{['Regular Spacing Gross',newline,'R^2 = ',num2str(round(gofReg210gross.rsquare,3))],...
                    ['Regular Spacing Net',newline,'R^2 = ',num2str(round(gofReg210net.rsquare,3))],...
                    ['sandTES Spacing Gross',newline,'R^2 = ',num2str(round(gofST210gross.rsquare,3))],...
                    ['sandTES Spacing Net',newline,'R^2 = ',num2str(round(gofST210net.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YLim=[100,350];
ax.YGrid='on';
saveas(fig,'Figure10.tiff');

save('resultsMICRO.mat','tab210','fitReg210gross','fitST210gross','wmf210','-append');


%% Comparison between 87 µm and 210 µm
%Create figure
fig=figure(11);
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
saveas(fig,'Figure11.tiff');


%% Analyze accuracies
tab87.P_FEED=[];
tab=[tab87;tab210];
names={'87mum_RegularPitch','87mum_sandTESPitch','210mum_RegularPitch','210mum_sandTESPitch'};
longNames={'d_p=87 µm, Regular Pitch','d_p=87 µm, sandTES Pitch','d_p=210 µm, Regular Pitch','d_p=210 µm, sandTES Pitch'};


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

tab.alpha_grossLower=tab.alpha_grossLower./tab.alpha_gross;
tab.alpha_grossUpper=tab.alpha_grossUpper./tab.alpha_gross;

tab.alpha_netLower=tab.alpha_netLower./tab.alpha_net;
tab.alpha_netUpper=tab.alpha_netUpper./tab.alpha_net;


%Write table and create figures
for i=1:length(names)
    ind=contains(tab.Dataset,names{i});
    
    tabAcc.FGmin(i)=min([1-tab.FGlower(ind),tab.FGupper(ind)-1],[],'all')*100;
    tabAcc.alpha_grossMin(i)=min([1-tab.alpha_grossLower(ind),tab.alpha_grossUpper(ind)-1],[],'all')*100;
    tabAcc.alpha_netMin(i)=min([1-tab.alpha_netLower(ind),tab.alpha_netUpper(ind)-1],[],'all')*100;
    
    tabAcc.FGmean(i)=mean([1-tab.FGlower(ind),tab.FGupper(ind)-1],'all')*100;
    tabAcc.alpha_grossMean(i)=mean([1-tab.alpha_grossLower(ind),tab.alpha_grossUpper(ind)-1],'all')*100;
    tabAcc.alpha_netMean(i)=mean([1-tab.alpha_netLower(ind),tab.alpha_netUpper(ind)-1],'all')*100;
    
    tabAcc.FGmax(i)=max([1-tab.FGlower(ind),tab.FGupper(ind)-1],[],'all')*100;
    tabAcc.alpha_grossMax(i)=max([1-tab.alpha_grossLower(ind),tab.alpha_grossUpper(ind)-1],[],'all')*100;
    tabAcc.alpha_netMax(i)=max([1-tab.alpha_netLower(ind),tab.alpha_netUpper(ind)-1],[],'all')*100;
    
    
    %Create figure
    fig=figure(100+i);
    clf(fig);
    ax=gca;
    colors=ax.ColorOrder;
    xlim=[1,nnz(ind)];
    hold on
    
    plot([tab.FGlower(ind),tab.FGupper(ind)],'Color',colors(1,:),'LineStyle','-');
    plot([tab.alpha_grossLower(ind),tab.alpha_grossUpper(ind)],'Color',colors(2,:),'LineStyle','--');
    plot([tab.alpha_netLower(ind),tab.alpha_netUpper(ind)],'Color',colors(3,:),'LineStyle',':','LineWidth',1);
    
    
    %Figure formatting
    legItems=repmat(line(),3,1);
    legItems(1)=plot(NaN,'Color',colors(1,:),'LineStyle','-');
    legItems(2)=plot(NaN,'Color',colors(2,:),'LineStyle','--');
    legItems(3)=plot(NaN,'Color',colors(3,:),'LineStyle',':','LineWidth',1);
    hold off
    
    title(ax,['Measurement Uncertainties, MICRO, ',longNames{i}]);
    xlabel(ax,'Dataset Index (-)');
    ylabel(ax,'Relative Uncertainty (-)');
    legend(legItems,{'FG','\alpha_{gross}','\alpha_{net}'},'Location','bestoutside');
    
    fig.Units='centimeters';
    fig.Position=[10,5,17,8.5];
    ax.XLim=xlim;
    ax.YGrid='on';
    saveas(fig,['Figure',num2str(100+i),'.tiff']);  
end





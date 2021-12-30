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
%https://doi.org/10.5281/zenodo.5802409
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5802407
%
%
%
%This script analyzes the LINI test rig data and creates all published
%figures.
%
%
%Requires all LINI data files in the same folder and all auxiliary
%functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
%   - accLINI.m
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

DeltaH_eps=50e-3;

d_p=146e-6;
rho_p=2650;
rho_Anorm=1.293;
g=9.81;
p_amb=1013.25e2;


dirCont=dir();


%% Plain tubes
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'LINI_Plain_'));


%Initialize table
varnames={'Dataset','I','TAin','TAout','T_bed',...
            'T_surf','U','VdotRota',...
            'dpFloor','pAin','pRota','p_eps1','p_eps2',...
            'phi','Vdot','eps'...
            'FG','FGlower','FGupper',...
            'alpha_gross','alpha_grossLower','alpha_grossUpper',...
            'alpha_net','alpha_netLower','alpha_netUpper'};
tabPlain=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tabPlain.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tabPlain{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,3:6}=tabloc{:,3:6}+273.15;
    tabloc{:,[8,15]}=tabloc{:,[8,15]}./60^2;
    tabloc.phi=deg2rad(tabloc.phi);
    
    tabPlain{i,[2:11,14:15]}=mean(tabloc{:,[2:11,14:15]});
    
    excludeP=tabloc.p_eps1<=0 | tabloc.p_eps1>=1000 | tabloc.p_eps2<=0 | tabloc.p_eps2>=1000;
    tabPlain{i,12:13}=mean(tabloc{~excludeP,12:13});
end


%Calculate bed porosities
eps1=1-tabPlain.p_eps1./(rho_p*g*DeltaH_eps);
eps2=1-tabPlain.p_eps2./(rho_p*g*DeltaH_eps);
tabPlain.eps=mean([eps1,eps2],2);


%Calculate degree of fluidization
p_Aplain=p_amb+tabPlain.pAin-tabPlain.dpFloor+rho_p*g*35e-3*(1-tabPlain.eps);
T_Aplain=mean([tabPlain.TAin,tabPlain.TAout],2);
rho_A=DryAir.rho(p_Aplain,T_Aplain);
Vdot=tabPlain.Vdot.*rho_Anorm./rho_A;
wmfPlain=w_mf(d_p,rho_p,p_Aplain,T_Aplain);
Vdot_mf=wmfPlain.*A_FB;
tabPlain.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
tabPlain.alpha_gross=tabPlain.U.*tabPlain.I.*cos(tabPlain.phi)./A./(tabPlain.T_surf-tabPlain.T_bed);


%Calculate net heat transfer coefficient
mDotA_loss=Vdot.*rho_A.*x_heated;
QdotLoss=mDotA_loss.*(DryAir.h(tabPlain.TAout)-DryAir.h(tabPlain.TAin));
tabPlain.alpha_net=(tabPlain.U.*tabPlain.I.*cos(tabPlain.phi)-QdotLoss)./A./(tabPlain.T_surf-tabPlain.T_bed);


%Calculate accuracy
tabPlain=accLINI(tabPlain);


%Create figure
fig=figure(12);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,5.5];
x=tabPlain.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tabPlain.alpha_gross;
scatter(ax,x,y,10,colors(1,:),'o');

[fitPlainGross,gofPlainGross]=fit(x,y,'poly2');
plot(xfit,fitPlainGross(xfit),'Color',colors(1,:),'LineStyle','-');


y=tabPlain.alpha_net;
scatter(ax,x,y,10,colors(2,:),'+');

[fitPlainNet,gofPlainNet]=fit(x,y,'poly2');
plot(xfit,fitPlainNet(xfit),'Color',colors(2,:),'LineStyle','--');


%Figure formatting
legItems=repmat(line(),2,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
hold off

title(ax,['HTCs LINI Rig, Plain Tube, w_{mf}=',num2str(round(mean(wmfPlain*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(legItems,{['Gross, R^2 = ',num2str(round(gofPlainGross.rsquare,3))],...
                    ['Net, R^2 = ',num2str(round(gofPlainNet.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure12.tiff');

save('resultsLINI.mat','tabPlain','fitPlainGross','p_Aplain','T_Aplain','wmfPlain');


%% Finned tubes
%Retrieve filenames
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'LINI_Finned_'));


%Initialize table
tabFinned=table('Size',[length(files),length(varnames)],'VariableTypes',[{'string'},repmat({'double'},1,length(varnames)-1)]);
tabFinned.Properties.VariableNames=varnames;


%Read files, transform to SI base units and calculate mean values
for i=1:length(files)
    tabFinned{i,1}={files{i}(1:end-4)};
    
    tabloc=readtable(files{i},'VariableNamingRule','preserve');
    tabloc{:,3:6}=tabloc{:,3:6}+273.15;
    tabloc{:,[8,15]}=tabloc{:,[8,15]}./60^2;
    tabloc.phi=deg2rad(tabloc.phi);
    
    tabFinned{i,[2:11,14:15]}=mean(tabloc{:,[2:11,14:15]});
    
    excludeP=tabloc.p_eps1<=0 | tabloc.p_eps1>=1000 | tabloc.p_eps2<=0 | tabloc.p_eps2>=1000;
    if nnz(excludeP)>height(tabloc)/2
        error('nnz(excludeP)>height(tabloc)/2');
    end
    tabFinned{i,12:13}=mean(tabloc{~excludeP,12:13});
end


%Calculate bed porosities
eps1=1-tabFinned.p_eps1./(rho_p*g*DeltaH_eps);
eps2=1-tabFinned.p_eps2./(rho_p*g*DeltaH_eps);
tabFinned.eps=mean([eps1,eps2],2);


%Calculate degree of fluidization
p_Afinned=p_amb+tabFinned.pAin-tabFinned.dpFloor+rho_p*g*35e-3*(1-tabFinned.eps);
T_Afinned=mean([tabFinned.TAin,tabFinned.TAout],2);
rho_A=DryAir.rho(p_Afinned,T_Afinned);
Vdot=tabFinned.Vdot.*rho_Anorm./rho_A;
wmfFinned=w_mf(d_p,rho_p,p_Afinned,T_Afinned);
Vdot_mf=wmfFinned.*A_FB;
tabFinned.FG=Vdot./Vdot_mf;


%Calculate gross heat transfer coefficient
tabFinned.alpha_gross=tabFinned.U.*tabFinned.I.*cos(tabFinned.phi)./A./(tabFinned.T_surf-tabFinned.T_bed);


%Calculate net heat transfer coefficient
mDotA_loss=Vdot.*rho_A.*x_heated;
QdotLoss=mDotA_loss.*(DryAir.h(tabFinned.TAout)-DryAir.h(tabFinned.TAin));
tabFinned.alpha_net=(tabFinned.U.*tabFinned.I.*cos(tabFinned.phi)-QdotLoss)./A./(tabFinned.T_surf-tabFinned.T_bed);


%Calculate accuracy
tabFinned=accLINI(tabFinned);


%Create figure
fig=figure(13);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,5.5];
x=tabFinned.FG;
xfit=linspace(xlim(1),xlim(2),50);
hold on


y=tabFinned.alpha_gross;
scatter(ax,x,y,10,colors(1,:),'o');

[fitFinnedGross,gofFinnedGross]=fit(x,y,'poly2');
plot(xfit,fitFinnedGross(xfit),'Color',colors(1,:),'LineStyle','-');


y=tabFinned.alpha_net;
scatter(ax,x,y,10,colors(2,:),'+');

[fitFinnedNet,gofFinnedNet]=fit(x,y,'poly2');
plot(xfit,fitFinnedNet(xfit),'Color',colors(2,:),'LineStyle','--');


%Figure formatting
legItems=repmat(line(),2,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:),'LineStyle','-','Marker','o');
legItems(2)=plot(NaN,NaN,'Color',colors(2,:),'LineStyle','--','Marker','+');
hold off

title(ax,['Virtual HTCs LINI Rig, Finned Tube, w_{mf}=',num2str(round(mean(wmfFinned*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Virtual HTC (W/m²K)');
legend(legItems,{['Gross, R^2 = ',num2str(round(gofFinnedGross.rsquare,3))],...
                    ['Net, R^2 = ',num2str(round(gofFinnedNet.rsquare,3))]},...
                    'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=xlim;
ax.YGrid='on';
saveas(fig,'Figure13.tiff');

save('resultsLINI.mat','tabFinned','fitFinnedGross','wmfFinned','p_Afinned','T_Afinned','-append');


%% Comparison between plain and finned tubes
%Create figure
fig=figure(14);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
xlim=[1.5,5.5];
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
saveas(fig,'Figure14.tiff');


%% Analyze accuracies
tab=[tabPlain;tabFinned];
names={'LINI_Plain','LINI_Finned'};
longNames={'Plain Tubes','Finned Tubes'};


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
    fig=figure(200+i);
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
    
    title(ax,['Measurement Uncertainties, LINI, ',longNames{i}]);
    xlabel(ax,'Dataset Index (-)');
    ylabel(ax,'Relative Uncertainty (-)');
    legend(legItems,{'FG','\alpha_{gross}','\alpha_{net}'},'Location','bestoutside');
    
    fig.Units='centimeters';
    fig.Position=[10,5,17,8.5];
    ax.XLim=xlim;
    ax.YGrid='on';
    saveas(fig,['Figure',num2str(200+i),'.tiff']);  
end




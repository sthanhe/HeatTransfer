%% Discussion results
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
%This script analyzes the results of the paper and creates all figures
%published in the "Discussion" section.
%
%
%Requires all auxiliary functions on the MATLAB path. All data files can be
%created with the respective "Analyze_..." scripts for the MICRO, LINI and
%TRINI experiments
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary classes and functions:
%   - @DryAir
%   - @SiO2
%   - andeenGlicksman.m
%   - epsilon.m
%   - gelperinEinstein.m
%   - grewal.m
%   - molerus.m
%   - zabrodsky.m
%   - w_mf.m
%Data files:
%   - resultsLINI.mat
%   - resultsMICRO.mat
%   - resultsTRINI.mat


%% General boundary conditions
d_tube=25e-3;

pitch_horzST=2.*d_tube;
pitch_vertST=2.5.*d_tube;
pitchST=mean([pitch_horzST,pitch_vertST]);

pitch_horzReg=3.1.*d_tube;
pitch_vertReg=3.1.*d_tube;
pitchReg=mean([pitch_horzReg,pitch_vertReg]);

d_p=146e-6;
d_p1=87.141e-6;
d_p2=210e-6;

eps_mf=0.45;    %bed porosity at minimum fluidization
eps_R=0.9;      %emissivity tube-bed
lambda_p=3;     %thermal conductivity of SiO2 
rho_p=2650;     %density of SiO2


%% Comparison of MICRO results with literature
load('resultsMICRO.mat');

%Boundary conditions
p_A=1013.25e2+4335;


%% MICRO 87 µm
%Correlations
FGlim=[3,7];
FG=linspace(FGlim(1),FGlim(2),height(tab87));

grewalST87=grewal(d_p1,rho_p,FG',d_tube,pitchST,p_A,tab87.T_BED1);
grewalReg87=grewal(d_p1,rho_p,FG',d_tube,pitchReg,p_A,tab87.T_BED1);
martin87=martin(d_p1,rho_p,lambda_p,FG',p_A,tab87.T_BED1,eps_mf,eps_R);
molerus87=molerus(rho_p,p_A,tab87.T_BED1);


%Create figure
fig=figure(23);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
hold on


%Plot fits
plot(FG,fitST87gross(FG),'Color',colors(1,:),'LineStyle','-');
plot(FG,fitReg87gross(FG),'Color',colors(1,:),'LineStyle','--');
plot(FG,grewalST87,'Color',colors(2,:),'LineStyle',':','LineWidth',1);
plot(FG,grewalReg87,'Color',colors(2,:),'LineStyle','-.');
plot(FG,martin87,'Color',colors(3,:),'LineStyle','-');
plot(FG,molerus87,'Color',colors(4,:),'LineStyle','--');
hold off


%Figure formatting
pos=round(length(FG)/4);
set(ax.Children(1),{'Marker','MarkerIndices','MarkerSize'},{'o',pos:pos:3*pos,3});
set(ax.Children(2),{'Marker','MarkerIndices','MarkerSize'},{'+',pos:pos:3*pos,3});

title(ax,['HTCs MICRO Rig, d_p=87 µm, w_{mf}=',num2str(round(mean(wmf87*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,{'sandTES Spacing','Regular Spacing','Grewal sandTES','Grewal Regular','Martin','Molerus'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=FGlim;
ax.YGrid='on';
saveas(fig,'Figure23.tiff');


%% MICRO 210 µm
%Correlations
FGlim=[1,6];
FG=linspace(FGlim(1),FGlim(2),height(tab210));

grewalST210=grewal(d_p2,rho_p,FG',d_tube,pitchST,p_A,tab210.T_BED1);
grewalReg210=grewal(d_p2,rho_p,FG',d_tube,pitchReg,p_A,tab210.T_BED1);
martin210=martin(d_p2,rho_p,lambda_p,FG',p_A,tab210.T_BED1,eps_mf,eps_R);
molerus210=molerus(rho_p,p_A,tab210.T_BED1);


%Create figure
fig=figure(24);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
hold on


%Plot fits
plot(FG,fitST210gross(FG),'Color',colors(1,:),'LineStyle','-');
plot(FG,fitReg210gross(FG),'Color',colors(1,:),'LineStyle','--');
plot(FG,grewalST210,'Color',colors(2,:),'LineStyle',':','LineWidth',1);
plot(FG,grewalReg210,'Color',colors(2,:),'LineStyle','-.');
plot(FG,martin210,'Color',colors(3,:),'LineStyle','-');
plot(FG,molerus210,'Color',colors(4,:),'LineStyle','--');
hold off


%Figure formatting
pos=round(length(FG)/4);
set(ax.Children(1),{'Marker','MarkerIndices','MarkerSize'},{'o',pos:pos:3*pos,3});
set(ax.Children(2),{'Marker','MarkerIndices','MarkerSize'},{'+',pos:pos:3*pos,3});

title(ax,['HTCs MICRO Rig, d_p=210 µm, w_{mf}=',num2str(round(mean(wmf210*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,{'sandTES Spacing','Regular Spacing','Grewal sandTES','Grewal Regular','Martin','Molerus'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=FGlim;
ax.YGrid='on';
saveas(fig,'Figure24.tiff');


%% Comparison of LINI results with literature
load('resultsLINI.mat');


%Correlations
FGlim=[1.5,5.5];
FG=linspace(FGlim(1),FGlim(2),height(tabPlain));

grewalLINI=grewal(d_p,rho_p,FG',d_tube,pitchST,p_Aplain,T_Aplain);
martinLINI=martin(d_p,rho_p,lambda_p,FG',p_Aplain,T_Aplain,eps_mf,eps_R);
molerusLINI=molerus(rho_p,p_Aplain,T_Aplain);


%Create figure
fig=figure(25);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
hold on


%Plot fits
plot(FG,fitPlainGross(FG),'Color',colors(1,:),'LineStyle','-');
plot(FG,grewalLINI,'Color',colors(2,:),'LineStyle','--');
plot(FG,martinLINI,'Color',colors(3,:),'LineStyle',':','LineWidth',1);
plot(FG,molerusLINI,'Color',colors(4,:),'LineStyle','-.');
hold off


title(ax,['(Virtual) HTCs LINI Rig, w_{mf}=',num2str(round(mean([wmfFinned;wmfPlain]*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,{'LINI','Grewal','Martin','Molerus'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=FGlim;
ax.YGrid='on';
saveas(fig,'Figure25.tiff');


%% Comparison of TRINI results with literature
load('resultsTRINI.mat');


%Correlations
FGlim=[1.5,4.5];
FG=linspace(FGlim(1),FGlim(2),height(tabPlainWithout)+height(tabPlainWith));

grewalTRINI=grewal(d_p,rho_p,FG',d_tube,pitchST,[p_AplainWithout;p_AplainWith],[T_AplainWithout;T_AplainWith]);
martinTRINI=martin(d_p,rho_p,lambda_p,FG',[p_AplainWithout;p_AplainWith],[T_AplainWithout;T_AplainWith],eps_mf,eps_R);
molerusTRINI=molerus(rho_p,[p_AplainWithout;p_AplainWith],[T_AplainWithout;T_AplainWith]);


%Create figure
fig=figure(26);
clf(fig);
ax=gca;
colors=ax.ColorOrder;
hold on


%Plot fits
plot(FG,fitGrossWithout(FG),'Color',colors(1,:),'LineStyle','-');
plot(FG,fitGrossWith(FG),'Color',colors(1,:),'LineStyle','--');
plot(FG,grewalTRINI,'Color',colors(2,:),'LineStyle',':','LineWidth',1);
plot(FG,martinTRINI,'Color',colors(3,:),'LineStyle','-.');
plot(FG,molerusTRINI,'Color',colors(4,:),'LineStyle','-');
hold off


%Figure formatting
pos=round(length(FG)/4);
set(ax.Children(1),{'Marker','MarkerIndices','MarkerSize'},{'o',pos:pos:3*pos,3});

title(ax,['HTCs TRINI Rig, w_{mf}=',num2str(round(mean([wmfplainWithout;wmfplainWith]*10^3),1)),' mm/s']);
xlabel(ax,'Fluidization Degree (-)');
ylabel(ax,'Heat Transfer Coefficient (W/m²K)');
legend(ax,{'Without Baffle','With Baffle','Grewal','Martin','Molerus'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
ax.XLim=FGlim;
ax.YGrid='on';
saveas(fig,'Figure26.tiff');


%% Air mass flow fluidization requirement
%Boundary conditions
T0=40+273.15;
T_A=T0:0.1:550+273.15;
p_A=1e5;

rho_A0=DryAir.rho(p_A,T_A(1));
wmf0=w_mf(d_p,rho_p,p_A,T_A(1));


%Relative required massflow for fluidization when using 146 µm particles,
%compared to 40°C
rho_A=DryAir.rho(p_A,T_A);
wmf=w_mf(d_p,rho_p,p_A,T_A);
mDotRel=wmf.*rho_A./(wmf0.*rho_A0);


%Create figure
fig=figure(27);
clf(fig);
ax=gca;

plot(ax,T_A-273.15,mDotRel,'LineStyle','-');


%Figure formatting
title(ax,'Air Mass Flow Fluidization Requirement');
xlabel(ax,'Air or Bed Temperature (°C)');
ylabel(ax,'$$\frac{\dot m_A (T)}{\dot m_A (40^\circ C)}$$','Interpreter','latex');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
xlim(ax,[T_A(1),T_A(end)]-273.15);
ax.YGrid='on';
saveas(fig,'Figure27.tiff');


%% Temperature dependence of different heat transfer correlations
%Boundary conditions
FG=4;
p_A=1e5;


%Temperature range
T_A=40:0.01:550;
T_A=T_A+273.15;


%Calculate heat transfer coefficients from correlations
alpha=NaN(length(T_A),6);
alpha(:,1)=andeenGlicksman(d_p,rho_p,FG,d_tube,p_A,T_A);
alpha(:,2)=grewal(d_p,rho_p,FG,d_tube,pitchST,p_A,T_A);
alpha(:,3)=molerus(rho_p,p_A,T_A);
alpha(:,4)=zabrodsky(d_p,rho_p,T_A);
alpha(:,5)=martin(d_p,rho_p,lambda_p,FG,p_A,T_A,eps_mf,eps_R);
alpha(:,6)=gelperinEinstein(d_p,rho_p,d_tube,pitch_horzST,pitch_vertST,p_A,T_A);


%Normalize to 40°C
alphaRel=alpha./alpha(1,:);


%Create figure
fig=figure(28);
clf(fig)
ax=gca;

plot(T_A-273.15,alphaRel);


%Figure formatting
set(ax.Children,{'LineStyle'},fliplr({'-';'--';':';'-.';'-';'--'}));
ax.Children(3).LineWidth=1;
pos=round(length(T_A)/4);
set(ax.Children(1),{'Marker','MarkerIndices','MarkerSize'},{'o',pos:pos:3*pos,3});
set(ax.Children(2),{'Marker','MarkerIndices','MarkerSize'},{'+',pos:pos:3*pos,3});


title(ax,'HTC over Temperature');
xlabel(ax,'Temperature (°C)');
ylabel(ax,'$$\frac{HTC (T)}{HTC (40^\circ C)}$$','Interpreter','latex');
legend(ax,{'Andeen / Glicksman','Grewal','Molerus',...
            'Zabrodsky','Martin','Gelperin / Einstein'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
xlim(ax,[T_A(1),T_A(end)]-273.15);
ax.YGrid='on';
saveas(fig,'Figure28.tiff');





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
%https://doi.org/10.5281/zenodo.5474020
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.5500329
%
%
%This script analyzes the results of the paper and creates all figures
%publishe in the "Discussion" section
%Requires all auxiliary functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - @DryAir
%   - @SiO2
%   - andeenGlicksman.m
%   - epsilon.m
%   - gelperinEinstein.m
%   - grewal.m
%   - molerus.m
%   - zabrodsky.m
%   - w_mf.m


%% Air mass flow fluidization requirement
%Temperature range
T=40:0.1:550;
T=T+273.15;


%Relative required massflow for fluidization when using 146 µm particles,
%compared to 40°C
y=w_mf(146e-6,2650,T).*DryAir.rho(T)./(w_mf(146e-6,2650,T(1)).*DryAir.rho(T(1)));


%Create figure
fig=figure(22);
clf(fig);
ax=gca;

plot(ax,T-273.15,y,'LineStyle','-');


%Figure formatting
title(ax,'Air Mass Flow Fluidization Requirement');
xlabel(ax,'Air or Bed Temperature (°C)');
ylabel(ax,'$$\frac{\dot m_A (T)}{\dot m_A (40^\circ C)}$$','Interpreter','latex');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
xlim(ax,[T(1),T(end)]-273.15);
ax.YGrid='on';
saveas(fig,'Figure22.tiff');


%% Temperature dependency of different heat transfer correlations
%Temperature range
T=40:0.01:550;
T=T+273.15;


%Geometry boundary conditions
d_tube=25e-3;                           %Tube diameter
pitch_horz=2.*d_tube;                   %Horizontal pitch
pitch_vert=2.5.*d_tube;                 %Vertical pitch
pitch=mean([pitch_horz,pitch_vert]);    %Mean pitch
FG=4;                                   %Degree of fluidization

d_p=146e-6;     %Particle diameter
rho_p=2650;     %Particle density (SiO2)


%Calculate heat transfer coefficients from correlations
alpha=NaN(length(T),5);
alpha(:,1)=andeenGlicksman(d_p,rho_p,FG,d_tube,T);
alpha(:,2)=grewal(d_p,rho_p,FG,d_tube,pitch,T);
alpha(:,3)=molerus(rho_p,T);
alpha(:,4)=zabrodsky(d_p,rho_p,T);
alpha(:,5)=gelperinEinstein(d_p,rho_p,d_tube,pitch_horz,pitch_vert,T);


%Normalize to 40°C
alphanorm=alpha./alpha(1,:);


%Create figure
fig=figure(23);
clf(fig)
ax=gca;

plot(T-273.15,alphanorm);


%Figure formatting
set(ax.Children,{'LineStyle'},{'-';'--';':';'-.';'-'});
ax.Children(3).LineWidth=1;
pos=round(length(T)/4);
set(ax.Children(1),{'Marker','MarkerIndices','MarkerSize'},{'o',pos:pos:3*pos,3});


title(ax,'HTC over Temperature');
xlabel(ax,'Temperature (°C)');
ylabel(ax,'$$\frac{HTC (T)}{HTC (40^\circ C)}$$','Interpreter','latex');
legend(ax,{'Andeen / Glicksman','Grewal','Molerus',...
            'Zabrodsky','Gelperin / Einstein'},'Location','bestoutside');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];
xlim(ax,[T(1),T(end)]-273.15);
ax.YGrid='on';
saveas(fig,'Figure23.tiff');





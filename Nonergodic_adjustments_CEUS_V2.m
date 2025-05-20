clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%##########################################################################
% Nonergodic adjustments for CEUS
%##########################################################################
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculate PSA (g) (ergodic and Nonergodic) for selected,
% Event location,Station location, Magnitude, Distance, and spectral Period based on
% nonergodic adjustments by Davatgari-Tafreshi and Pezeshk (2025):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Please email us if you find any errors in the code or have any suggestions. 
% Thank you!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set desired M, R, T arrays
%   txt file of the GMM coefficients put in INPUT directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ergodic and nonergodic PSA (g) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-Loading the Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_m=readtable('Mean_Site.csv',"VariableNamingRule","preserve"); %Site adjustments ds2s
site_e=readtable('Epistemic_Site.csv',"VariableNamingRule","preserve"); %Epistemic uncertainty of Site adjustments SE(ds2s)
source_m=readtable('Source_Mean.csv',"VariableNamingRule","preserve");%Source adjustments dl2l
source_e=readtable('Source_Epistemic.csv',"VariableNamingRule","preserve");%Epistemic uncertainty of Source adjustments SE(dl2l)
path_m=readtable('Path_Mean.csv',"VariableNamingRule","preserve");%Path adjustments dp2p
path_e=readtable('Path_Epistemic2.csv',"VariableNamingRule","preserve");%Epistemic uncertainty of Source adjustments SE(dp2p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_rup_comb=50;%Distance
Mw_comb=7;%Magnitude
%For stations that are not on the list.
CP_inout='outside';%'inside'
Sed_depth=10000;
size1=20;% plot font size

Station_lat=34.5454;%Station latitude   %Station_lat = str2double(cell2mat(inputdlg('Please enter the latitude of the station from the station list:', 'Value Input', [1 50])));
Station_lon=-93.5765;%Station longitude %Station_lon = str2double(cell2mat(inputdlg('Please enter the longitude of the station from the station list:', 'Value Input', [1 50])));
Source_lat =37;%Source cell latitude	%Source_lat = str2double(cell2mat(inputdlg('Please enter the latitude of the Source from the station list:', 'Value Input', [1 50])));
Source_lon=-98;%Source cell longitude   %Source_lon = str2double(cell2mat(inputdlg('Please enter the longitude of the Source from the station list:', 'Value Input', [1 50])));

Period = [ 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10];
tau_l2l=[0.435000000000000	0.446000000000000	0.464000000000000	0.473000000000000	0.484000000000000	0.489000000000000	0.502000000000000	0.473000000000000	0.434000000000000	0.421000000000000	0.418000000000000	0.414000000000000	0.428000000000000	0.431000000000000	0.444000000000000	0.451000000000000	0.483000000000000	0.496000000000000	0.468000000000000	0.450000000000000	0.433000000000000	0.412000000000000	0.416000000000000];
phi_p2p=[0.423000000000000	0.632000000000000	0.585000000000000	0.562000000000000	0.558000000000000	0.555000000000000	0.515000000000000	0.476000000000000	0.405000000000000	0.361000000000000	0.345000000000000	0.345000000000000	0.357000000000000	0.360000000000000	0.355000000000000	0.348000000000000	0.346000000000000	0.346000000000000	0.369000000000000	0.384000000000000	0.395000000000000	0.437000000000000	0.484000000000000];
phi_s2s_out=[0.697000000000000	0.780000000000000	0.779000000000000	0.772000000000000	0.772000000000000	0.771000000000000	0.752000000000000	0.755000000000000	0.719000000000000	0.696000000000000	0.691000000000000	0.682000000000000	0.653000000000000	0.633000000000000	0.579000000000000	0.541000000000000	0.512000000000000	0.469000000000000	0.467000000000000	0.475000000000000	0.479000000000000	0.486000000000000	0.501000000000000];
phi_s2s_in=[0.443000000000000	0.415000000000000	0.413000000000000	0.414000000000000	0.422000000000000	0.428000000000000	0.445000000000000	0.455000000000000	0.463000000000000	0.471000000000000	0.472000000000000	0.499000000000000	0.511000000000000	0.526000000000000	0.506000000000000	0.503000000000000	0.511000000000000	0.464000000000000	0.493000000000000	0.499000000000000	0.468000000000000	0.432000000000000	0.430000000000000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3-Find earthquakes cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon_lat_polygon = load('lon_lat_cp.txt');
% Source location
A = [Source_lat, Source_lon]; % [longitude, latitude]
% Define the latitude and longitude ranges
lat_min = 24; % Minimum latitude
lat_max =  53; % Maximum latitude
lon_min = -107; % Minimum longitude
lon_max =  -58; % Maximum longitude
% Define the grid resolution
grid_resolution = 2.0;
% Grid edges
lat_edges = lat_min:grid_resolution:lat_max;
lon_edges = lon_min:grid_resolution:lon_max;
% Find the index of the lat and lon in their respective ranges
lat_idx = find(lat_edges <= Source_lat, 1, 'last');
lon_idx = find(lon_edges <= Source_lon, 1, 'last');
% Check if the found indices are within bounds
if isempty(lat_idx) || isempty(lon_idx)
    error('Source location is out of defined region bounds.');
end
% Determine the subregion number
subregion_idx = (lat_idx - 1) * length(lon_edges) + lon_idx;
% Display the subregion index
disp(['Source location is in subregion: ', num2str(subregion_idx)]);

% Create the latitude and longitude vectors with the specified resolution
lat_edges  = lat_min:grid_resolution:lat_max;
lon_edges = lon_min:grid_resolution:lon_max;

% Calculate the midpoints of the grid cells
lat_midpoints = lat_edges(1:end-1) + grid_resolution / 2;
lon_midpoints = lon_edges(1:end-1) + grid_resolution / 2;

%-Plot earthquakes cells
figure;
ax = axesm('MapProjection', 'mercator'); % Use Mercator projection or another appropriate projection
hold on;

LAT = lon_lat_polygon(:,2);
LON= lon_lat_polygon(:,1);
plotm(LAT, LON, 'k','LineWidth',2.5)
hold on
rp=[0.5 0.5 0.5];
for lon = lon_edges
    plotm([lat_min, lat_max], [lon, lon], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

end
for lat = lat_edges
    plotm([lat lat], [lon_min lon_max], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end
set(ax, 'FontSize', 12, 'FontName', 'Times');
hold on;
plotm(Source_lat, Source_lon, 'Marker', 'o', 'LineStyle', 'none', ...
    'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
hold on
plotm(lat_edges(lat_idx), lon_edges(lon_idx), 'Marker', 'o', 'LineStyle', 'none', ...
    'MarkerSize', 6, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
hold on
plotm(lat_midpoints(lat_idx), lon_midpoints(lon_idx), 'Marker', 'o', 'LineStyle', 'none', ...
    'MarkerSize', 6, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
% Load and plot state boundaries
states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(states, 'DisplayType', 'polygon', 'FaceColor', 'None', 'EdgeColor', 'k');
% Add labels or titles
xlabel('Longitude');
ylabel('Latitude');
axis tight;
xlim([-1.89 -1])

%Source Input
lat_source=source_m.Latitude;
lon_source=source_m.Longitude;

PGA_SA_R1=source_m(:,3:25);
PGA_SA_F1=table2array(PGA_SA_R1);
Source_Mean=[PGA_SA_F1];
clear PGA_SA_R1 PGA_SA_F1
PGA_SA_R1=source_e(:,3:25);
PGA_SA_F1=table2array(PGA_SA_R1);
Source_epi=[PGA_SA_F1];
clear PGA_SA_R1 PGA_SA_F1

%Source adjusment
for i=1:length(lat_source)
    if lat_source(i)==Source_lat && lon_source(i)==Source_lon

        hold on
        plotm(lat_source(i), lon_source(i), 'Marker', 'o', 'LineStyle', 'none', ...
            'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        Source_adj=Source_Mean(i,:);
        Source_adj_Epi=Source_epi(i,:);
    end
end

if isempty(Source_adj)==1
    Source_adj=zeros(1,23);
    Source_adj_Epi=tau_l2l;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4-Find station
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Site Input
lat_site=site_m.StationLatitude;
lon_site=site_m.StationLongitude;
Z_T=site_m.CoastalPlainsSedThick;
PGA_SA_R1=site_m(:,4:26);
PGA_SA_F1=table2array(PGA_SA_R1);
Site_Mean=[PGA_SA_F1];
clear PGA_SA_R1 PGA_SA_F1
PGA_SA_R1=site_e(:,4:26);
PGA_SA_F1=table2array(PGA_SA_R1);
Site_epi=[PGA_SA_F1];
clear PGA_SA_R1 PGA_SA_F1
%Site adjusment
for i=1:length(lat_site)
    if lat_site(i)==Station_lat && lon_site(i)==Station_lon

        hold on
        plotm(lat_site(i), lon_site(i), 'Marker', '^', 'LineStyle', 'none', ...
            'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
        Site_adj=Site_Mean(i,:);
        Site_adj_Epi=Site_epi(i,:);
        Z_Sed=Z_T(i);
    end
end
  plotm([lat_midpoints(lat_idx), Station_lat], ...
          [lon_midpoints(lon_idx), Station_lon], ...
          'Color', 'g', ...
          'LineWidth', 1.5);
print('-djpeg',['Fig2' '.jpg'],'-r300')

if isempty(Site_adj)==1
    Site_adj=zeros(1,23);

    if strcmp(CP_inout,'outside')==1
        Site_adj_Epi=phi_s2s_out;
        Z_Sed='NaN';
    else
        Site_adj_Epi=phi_s2s_in;
        Z_Sed=Sed_depth;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5-Find Path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Site Input
lat_event_P=path_m.EarthquakeLatitude;
lon_event_P=path_m.EarthquakeLongitude;
lat_station_P=path_m.StationLatitude;
lon_station_P=path_m.StationLongitude;

PGA_SA_R1=path_m(:,9:31);
PGA_SA_F1=table2array(PGA_SA_R1);
Path_Mean=[PGA_SA_F1];
clear PGA_SA_R1 PGA_SA_F1
PGA_SA_R1=path_e(:,9:31);
PGA_SA_F1=table2array(PGA_SA_R1);
Path_epi=[PGA_SA_F1];
clear PGA_SA_R1 PGA_SA_F1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6-Path adjusment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Processing... please wait.');
u=1;
for i=1:length(lat_event_P)
    % Find the index of the lat and lon in their respective ranges
    lat_idx_EP = find(lat_edges <= lat_event_P(i), 1, 'last');
    lon_idx_EP = find(lon_edges <= lon_event_P(i), 1, 'last');
    if lat_idx_EP==(lat_idx) && lon_idx_EP==(lon_idx)

        Path_adj1(u,:)=Path_Mean(i,:);
        Path_adj1_slat(u,1)=lat_station_P(i);
        Path_adj1_slon(u,1)=lon_station_P(i);
        Path_adj1_Epi(u,:)=Path_epi(i,:);
        u=u+1;
    else
        
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Processing... please wait.');

u=1;
for i=1:length(Path_adj1_slat)

    if Path_adj1_slat(i)==Station_lat && Path_adj1_slon(i)==Station_lon
        Path_adj(u,:)=Path_adj1(i,:);
        Path_adj_Epi(u,:)=Path_adj1_Epi(i,:);
        u=u+1;
    end
end

Path_adj_f=mean(Path_adj);
Path_adj_Epi_f=mean(Path_adj_Epi);

if isempty(Path_adj_f)==1
    Path_adj_f=zeros(1,23);
    Path_adj_Epi_f=phi_p2p;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7-GMM selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vector=readtable("vector.csv");%Akhani et al. (2024) for areas inside the Coastal Plain (ADP24)

C=table2array(vector);

gmpe_coeff=load('Coefficients.txt');%Pezeshk et al. (2018) for areas outside the Coastal Plain (PZ18)
if isnan(Z_Sed)==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Median PZ18 and ADP24
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:length(Period)
        for i=1:length(R_rup_comb)
            %PZ18
            R=sqrt(gmpe_coeff(11,j)^2+R_rup_comb(i)^2);
            log10_Y(i,j)=gmpe_coeff(1,j)+gmpe_coeff(2,j)*Mw_comb+gmpe_coeff(3,j)*Mw_comb^2+(gmpe_coeff(4,j)+gmpe_coeff(5,j)*Mw_comb)*min(log10(R),log10(60))+...
                (gmpe_coeff(6,j)+gmpe_coeff(7,j)*Mw_comb)*max(min(log10(R/60),log10(120/60)),0)+(gmpe_coeff(8,j)+gmpe_coeff(9,j)*Mw_comb)*max(log10(R/120),0)+...
                +gmpe_coeff(10,j)*R;

            PSA_erg(i,j)=log(10^(log10_Y(i,j)));

        end
    end

else
    for j=1:length(Period)
        for i=1:length(R_rup_comb)
            %PZ18
            R=sqrt(gmpe_coeff(11,j)^2+R_rup_comb(i)^2);
            log10_Y(i,j)=gmpe_coeff(1,j)+gmpe_coeff(2,j)*Mw_comb(i)+gmpe_coeff(3,j)*Mw_comb(i)^2+(gmpe_coeff(4,j)+gmpe_coeff(5,j)*Mw_comb(i))*min(log10(R),log10(60))+...
                (gmpe_coeff(6,j)+gmpe_coeff(7,j)*Mw_comb(i))*max(min(log10(R/60),log10(120/60)),0)+(gmpe_coeff(8,j)+gmpe_coeff(9,j)*Mw_comb(i))*max(log10(R/120),0)+...
                +gmpe_coeff(10,j)*R;
            %ADP24
            C1=C(1,j);
            C2=C(2,j);
            C3=C(3,j);
            C4=C(4,j);
            C5=C(5,j);
            Fcp=C1*(Z_sed(i)/10)^C2+C3*log(R_rup_comb(i))+C4*R_rup_comb(i)+C5;
            PSA_erg(i,j)=log(10^(log10_Y(i,j)+Fcp));


        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NGMM1_M=exp(PSA_erg+Site_adj);
NGMM2_M=exp(PSA_erg+Site_adj+Source_adj);
NGMM3_M=exp(PSA_erg+Site_adj+Source_adj+Path_adj_f);

NGMM1_e_u=exp(PSA_erg+Site_adj+Site_adj_Epi);
NGMM2_e_u=exp(PSA_erg+Site_adj+Source_adj+Site_adj_Epi+Source_adj_Epi);
NGMM3_e_u=exp(PSA_erg+Site_adj+Source_adj+Path_adj_f+Site_adj_Epi+Source_adj_Epi+Path_adj_Epi_f);
NGMM1_e_l=exp(PSA_erg+Site_adj-Site_adj_Epi);
NGMM2_e_l=exp(PSA_erg+Site_adj+Source_adj-Site_adj_Epi-Source_adj_Epi);
NGMM3_e_l=exp(PSA_erg+Site_adj+Source_adj+Path_adj_f-Site_adj_Epi-Source_adj_Epi-Path_adj_Epi_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8-PSA Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0.1 0.1 0.5 0.9])

subplot(2,2,1)
set(gca,'FontSize',14,'LineWidth',1.3);
h=loglog(Period,exp(PSA_erg),'-ok','LineWidth',1.5);
hold on
h1=loglog(Period,NGMM1_M,'-sr','LineWidth',1.5);
hold on
h2=loglog(Period,NGMM1_e_u,'--r','LineWidth',1.5);
hold on
h3=loglog(Period,NGMM1_e_l,'--r','LineWidth',1.5);
hold on
ylabel('PSA (g)','fontname','times','fontsize',size1);
grid on
a=get(h,'Parent');set(a,'XTick',[0.01 0.1 1 3 10], ...
    'XTickLabel',[0.01 0.1 1 3 10]);
legend('Erg GMM','Erg GMM+\deltaS2S','\psi','fontsize',13,'Location','southwest','NumColumns', 1)
set(gca,'fontsize',15,'linewidth',1.0);
ylim([10^-4 10^0])
xlim([0.01 10])
ax = gca; % Get the current axes
ax.OuterPosition(3) = ax.OuterPosition(3) + 0.06;
subplot(2,2,2)
set(gca,'FontSize',14,'LineWidth',1.3);
h=loglog(Period,exp(PSA_erg),'-ok','LineWidth',1.5);
hold on
h1=loglog(Period,NGMM2_M,'-s','LineWidth',2.5,'Color',[0.3, 0.6, 1]);
hold on
h2=loglog(Period,NGMM2_e_u,'--','LineWidth',2.5,'Color',[0.3, 0.6, 1]);
hold on
h3=loglog(Period,NGMM2_e_l,'--','LineWidth',2.5,'Color',[0.3, 0.6, 1]);
hold on
grid on
a=get(h,'Parent');set(a,'XTick',[0.01 0.1 1 3 10], ...
    'XTickLabel',[0.01 0.1 1 3 10]);
legend('Erg GMM','Erg GMM+\deltaS2S+\deltaL2L','\psi','fontsize',13,'Location','southwest','NumColumns', 1)
set(gca,'fontsize',15,'linewidth',1.0);
ylim([10^-4 10^0])
xlim([0.01 10])
ax = gca; % Get the current axes
ax.OuterPosition(3) = ax.OuterPosition(3) + 0.06;
subplot(2,2,3)
set(gca,'FontSize',14,'LineWidth',1.3);
h=loglog(Period,exp(PSA_erg),'-ok','LineWidth',1.5);
hold on
h1=loglog(Period,NGMM3_M,'-sg','LineWidth',2);
hold on
h2=loglog(Period,NGMM3_e_u,'--g','LineWidth',2);
hold on
h3=loglog(Period,NGMM3_e_l,'--g','LineWidth',2);
hold on
ylabel('PSA (g)','fontname','times','fontsize',size1);
xlabel('Period (s)','fontname','times','fontsize',size1);
grid on
a=get(h,'Parent');set(a,'XTick',[0.01 0.1 1 3 10], ...
    'XTickLabel',[0.01 0.1 1 3 10]);
legend('Erg GMM','Erg GMM+\deltaS2S+\deltaL2L+\deltaP2P','\psi','fontsize',13,'Location','southwest','NumColumns',1)
set(gca,'fontsize',15,'linewidth',1.0);
ylim([10^-4 10^0])
xlim([0.01 10])
ax = gca; % Get the current axes
ax.OuterPosition(3) = ax.OuterPosition(3) + 0.06;
ax.OuterPosition(2) = ax.OuterPosition(2) + 0.06;
%%%%%%

subplot(2,2,4)
set(gca,'FontSize',14,'LineWidth',1.3);
h1=loglog(Period,exp(PSA_erg),'-ok','LineWidth',1.5);
hold on
h2=loglog(Period,NGMM1_M,'-sr','LineWidth',1.5);
hold on
h3=loglog(Period,NGMM2_M,'-s','LineWidth',2.5,'Color',[0.3, 0.6, 1]);
hold on
h4=loglog(Period,NGMM3_M,'-sg','LineWidth',2);
hold on
legend('Erg GMM','Erg GMM+\deltaS2S','Erg GMM+\deltaS2S+\deltaL2L',...
    'Erg GMM+\deltaS2S+\deltaL2L+\deltaP2P','fontsize',13,'Location','southwest','NumColumns',1)
xlabel('Period (s)','fontname','times','fontsize',size1);
grid on
a=get(h1,'Parent');set(a,'XTick',[0.01 0.1 1 3 10], ...
    'XTickLabel',[0.01 0.1 1 3 10]);
set(gca,'fontsize',15,'linewidth',1.0);
ylim([10^-4 10^0])
xlim([0.01 10])
ax = gca; % Get the current axes
ax.OuterPosition(3) = ax.OuterPosition(3) + 0.06;
    ax.OuterPosition(2) = ax.OuterPosition(2) + 0.06;
print(gcf,'PSA.jpg', '-dpng','-r300')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
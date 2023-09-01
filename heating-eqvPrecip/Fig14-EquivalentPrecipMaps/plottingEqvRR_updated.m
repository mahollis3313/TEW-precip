%Old method: TEW_integ_heating was in K/day, then divided by lat/lon
%counter, get clim mean. Then, conversion to mm/day by *40, then *365. Wayyyy too high.
%New method: TEW_integ_heating in K/day for each year. Take climatological mean for 20-yr period.
%Take to K/hr (/24), then *3 for 3-hourly accumulations (per year). Then convert to mm by factor of 40.

close all; clear;
load TEW_sumLH.mat %the average of each yrly 2D integration
%load TEW_avgLH2.mat
LHavg_clim=LH_accum_avg_clim./24.*3; %K/day to K/hr*3hrs
LHavg_clim=LHavg_clim.*(1/9.81);
%avg_eqvRR=(LHavg_clim.*10.*1000./250).*365; %multiply by 365 to get ann avg accum in mm
avg_eqvRR=LHavg_clim.*10.*1000./250;
%GET GLOBAL MEAN CONDITIONAL MEAN
glob_mean = num2str(round((nansum(nansum(avg_eqvRR))))./...
    (sum(sum(~isnan(avg_eqvRR)))),3);


lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)'; load coast;

intervals = [0, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0, 250.0, 500.0, 1000.0];
numIntervals = numel(intervals);
% Determine the index for each interval
indexData = zeros(size(avg_eqvRR));
for i = 2:numIntervals
    indexData(avg_eqvRR >= intervals(i-1) & avg_eqvRR < intervals(i)) = i-1;
end

% Define the RGB values for each interval
colors = [
    0, 0.369, 0.02;       % Green shades
    0.373, 0.6, 0.18;       % Green shades
    0.565, 0.749, 0.02;       % Green shades
    0.82, 0.831, 0.024;       % Greenish-yellow shade
    0.996, 1, 0.659;       % Yellow shades
    0.969, 0.769, 0.275;       % Yellow shades
    %1.0, 0.6, 0;       % Orange shades
    0.969, 0.576, 0.204;       % Orange shades
    0.878, 0.278, 0.125;       % Red shades
    0.6,0.016,0.016        % Red shades (or 1.0,0.0,0)
];

close all;
figure('Position', [100, 100, 800, 600]);  % Adjust the figure size and position
pcolor(lon_tdt,lat_tdt,indexData'); shading flat;
colormap(colors);
caxis([-0.5, numIntervals-0.5]); % Set color limits to cover the full range
cb=colorbar;
cb.Location = 'southoutside';  cb.TickLabels = intervals'; cb.Scale = 'log';
cb.Label.String = 'mm'; %day^{-1}
cb.TickLength = 0; %Hide the default colorbar tick marks
hold on;
plot(long,lat,'k','LineWidth',.75) %plot coastlines
axis tight
axis equal
axis([ -179.75 179.75 -36.75 36.75])
xticklabels({['150' char(176) 'W'],['100' char(176) 'W'],['50' char(176) 'W'],['0' char(176)],...
    ['50' char(176) 'E'],['100' char(176) 'E'],['150' char(176) 'E']})
yticklabels({['20' char(176) 'S'],['0' char(176)],['20' char(176) 'N']})
title(['                                Annual Average Equivalent Surface Rain (2001-2020)                         \mu = ',glob_mean])
print('-f1','/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/annacumTEST_cb3','-dpng')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear;
load TEW_sumLH.mat %the average of each yrly 2D integration
%load TEW_avgLH2.mat
LHavg_clim=LH_accum_avg_clim./24.*3; %K/day to K/hr*3hrs
LHavg_clim=LHavg_clim.*(1/9.81);
%avg_eqvRR=(LHavg_clim.*10.*1000./250).*365; %multiply by 365 to get ann avg accum in mm
avg_eqvRR=LHavg_clim.*10.*1000./250;
%GET GLOBAL MEAN CONDITIONAL MEAN
glob_mean = num2str(round((nansum(nansum(avg_eqvRR))))./...
    (sum(sum(~isnan(avg_eqvRR)))),3);


lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)'; load coast;

c1=[1 1 1]; 
c2=[.94 1 1];
c3=[.68 .85 .90];
c4=[.55 .75 .84];
c5=[.27 .51 .71];
c6=[0 .65 .58];
c7=[.19 .70 .10];
c8=[.60 .80 .20];
c9=[.82 .89 .19];
c10=[1 .80 .20];
c11=[1 .47 0];
c12=[1 .22 0];
c13=[.89 .26 .20];
c14=[.77 .12 .23];
c15=[.62 .11 .20];
mycolormap=rgbmap(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,15);

close all;
figure(1)
pcolor(lon_tdt,lat_tdt,avg_eqvRR');
colormap(mycolormap); shading flat;
caxis([0 1000])
hold on;
plot(long,lat,'k','LineWidth',.75) %plot coastlines
axis tight
axis equal
cb=colorbar('horizonta','position',[0.2 0.3 0.6 0.01]); %left,bott,W,H
cb.Label.String = 'mm';
axis([ -179.75 179.75 -36.75 36.75])
xticklabels({['150' char(176) 'W'],['100' char(176) 'W'],['50' char(176) 'W'],['0' char(176)],...
    ['50' char(176) 'E'],['100' char(176) 'E'],['150' char(176) 'E']})
yticklabels({['20' char(176) 'S'],['0' char(176)],['20' char(176) 'N']})
title(['                    Annual Average Equivalent Surface Rain (2001-2020)               \mu = ',glob_mean])
print('-f1','/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/MER_EqvPcp_850waves-cb2','-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% repeat but using the TOTPRECORR
clear;
load TEW_avgPCP_clim_accum.mat;
%Get global mean conditional mean
glob_mean = num2str(round((nansum(nansum(PCPavg_clim))))./...
    (sum(sum(~isnan(PCPavg_clim)))),3);

lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)'; load coast;

c1=[1 1 1]; 
c2=[.94 1 1];
c3=[.68 .85 .90];
c4=[.55 .75 .84];
c5=[.27 .51 .71];
c6=[0 .65 .58];
c7=[.19 .70 .10];
c8=[.60 .80 .20];
c9=[.82 .89 .19];
c10=[1 .80 .20];
c11=[1 .47 0];
c12=[1 .22 0];
c13=[.89 .26 .20];
c14=[.77 .12 .23];
c15=[.62 .11 .20];
mycolormap=rgbmap(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,14);

close all;
figure(1)
pcolor(lon_tdt,lat_tdt,PCPavg_clim');
colormap(mycolormap); shading flat;
caxis([0 1000])
hold on;
plot(long,lat,'k','LineWidth',.75) %plot coastlines
axis tight
axis equal
cb=colorbar('horizonta','position',[0.2 0.3 0.6 0.01]); %left,bott,W,H
cb.Label.String = 'mm';
axis([ -179.75 179.75 -36.75 36.75])
xticklabels({['150' char(176) 'W'],['100' char(176) 'W'],['50' char(176) 'W'],['0' char(176)],...
    ['50' char(176) 'E'],['100' char(176) 'E'],['150' char(176) 'E']})
yticklabels({['20' char(176) 'S'],['0' char(176)],['20' char(176) 'N']})
title(['                   Annual Average Accumulated MER Precip (2001-2020)               \mu = ',glob_mean])
print('-f1','/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/MER_PRECTOTCORR_850waves','-dpng')

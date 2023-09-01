
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMERG period: 2001-2020
%3-hourly NH/SH wave tracks (850 hPa)
% REPEAT USING 700 hPa 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear;
lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)';

dirname = '/home/c770m986/TEW_MER_tracking/FilteredNoTCUpdated/850hPa/';
%for y = 2001:2020
    y=2019; year = num2str(y);
    Ntrack_fname=strcat(dirname,'post_filter.',year,'.NH.NoTC.nc');
    Strack_fname=strcat(dirname,'post_filter.',year,'.SH.NoTC.nc');
    
    %Read in the tracking info
    N_time=ncread(Ntrack_fname,'time');
    N_lon=ncread(Ntrack_fname,'longitude');
    N_lat=ncread(Ntrack_fname,'latitude');
    S_time=ncread(Strack_fname,'time');
    S_lon=ncread(Strack_fname,'longitude');
    S_lat=ncread(Strack_fname,'latitude');

    %Combine N.H. & S.H. arrays/variables
    time=cat(1,N_time,S_time); 
    lon=cat(1,N_lon,S_lon); lon=wrapTo180(lon);  
    lat=cat(1,N_lat,S_lat);
    time=datetime(1979,12,1,0,0,0) + hours(time);
    time.Format='yyyyMMddHHmmss';  
    formatOut='yyyymmddHHMMss'; 
    time=datestr(time,formatOut);
    clear N_time N_lon N_lat S_time S_lon S_lat

    %Initialize 
    TEW_integ_heating=zeros(length(lon_tdt),length(lat_tdt)); %x,y or (576 x 161)
    counter_integ = zeros(length(lon_tdt),length(lat_tdt));
    for i = 1:length(time) 
        %First temporal matching
        t_w = time(i,:);
        f_w = t_w(1:8);
        f_list = strcat('/data/reanalysis/MERRA2/Daily/DTDT_TRMM_zlvls_v2/MERRA2_asm.tdt.3hr.',f_w,'.nc');
        dtdt_time = ncread(f_list,'time');
        yr = t_w(9:10);
        yr_out = str2double(yr);
        t_mm = yr_out*60; %convert from UTC hour to UTC minutes
        t_index = find(dtdt_time==t_mm);
        LH_Vavg = ncread(f_list,'LH_Vavg'); 

        %Next, spatial matching using a 500-km radius or 5-degree radius
        lat_window1=lat(i)-5; lat_window2=lat(i)+5;
        lon_window1=lon(i)-5; lon_window2=lon(i)+5; 

        lat_diff1 = abs(lat_tdt - lat_window1); 
        lat_min1 = min(lat_diff1);
        lat_index1 = find(lat_diff1==lat_min1);
        lat_diff2 = abs(lat_tdt - lat_window2); 
        lat_min2 = min(lat_diff2);
        lat_index2 = find(lat_diff2==lat_min2);

        lon_diff1 = abs(lon_tdt - lon_window1); 
        lon_min1 = min(lon_diff1);
        lon_index1 = find(lon_diff1==lon_min1);
        lon_diff2 = abs(lon_tdt - lon_window2); 
        lon_min2 = min(lon_diff2);
        lon_index2 = find(lon_diff2==lon_min2);
        
        %Assign dtdtmst values to corresponding temporary array within radius
        temp_heating = zeros(576,161);
        temp_heating(lon_index1:lon_index2,lat_index1:lat_index2) = ...
            LH_Vavg(lon_index1:lon_index2,lat_index1:lat_index2,t_index);
        temp_heating(isnan(temp_heating)) = 0; %replace all NaN values with zeros, while keeping all other values the same
        TEW_integ_heating = TEW_integ_heating + temp_heating; %accumlated heating (already in K/day)
        
        %Need a counter variable that will keep track of the lat/lons that are added in the integration
        %initialize 'counter_integ' before start of first loop
        temp_counter = zeros(576,161);
        temp_counter(lon_index1:lon_index2,lat_index1:lat_index2) = 1;
        counter_integ = counter_integ + temp_counter; %'accumulated' counts (576x161)
    end         


    var_name1 = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/700hPa_waves/TEW_accumavgLH_',year,'.mat');
    save(var_name1,'TEW_integ_heating')      %already in K/day

    var_name = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/700hPa_waves/TEW_latlon_counter_',year,'.mat');
    save(var_name,'counter_integ')     


%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clear; close all; 

folder1 = '/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/';
filelist1 = dir([folder1 'TEW_accumavgLH_*.mat']);
filelist2 = dir([folder1 'TEW_latlon_counter_*.mat']);
LH_ann = zeros(576,161,20);                             
counter_ann = zeros(576,161,20); 
counter_ann_ALL = zeros(576,161,20); 
for i = 1:20
    load(strcat(filelist1(i).name));
    LH_ann(:,:,i) = TEW_integ_heating; %K/day
    load(strcat(filelist2(i).name));
    counter_ann(:,:,i) = counter_integ; %counter 
end

LH_accum_avg_clim = squeeze(nanmean(LH_ann,3)); %using nanmean or nansum gives same answer
%LH_accum_avg_clim2 = squeeze(nansum(LH_ann,3));
%TEW_counter_accum = squeeze(nansum(counter_ann,3));%also, using nanmean or mean gives same answer
counter_clim = squeeze(nanmean(counter_ann,3));

LHavg_clim = LH_accum_avg_clim ./ counter_clim;


v_name = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/TEW_counter_clim.mat');
save(v_name,'counter_clim')

var_name = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/TEW_avgLH.mat');
save(var_name,'LHavg_clim')
% 
% var_namev2=strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/TEW_avgLH2.mat');
% save(var_namev2,'LHavg_clim2')
% 
% var_namev3=strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/TEW_sumLH.mat');
% save(var_namev3,'LH_accum_avg_clim')

% var_namev4=strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/TEW_sumLH2.mat');
% save(var_namev4,'LH_accum_avg_clim2')

% var_name1 = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/IMERG_period/TEW_avgLH.mat');
% save(var_name1,'LHavg_clim')
% 
% var_name2 = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/IMERG_period/TEW_counter_clim2.mat');
% save(var_name2,'counter_clim')
% 
% var_name3 = strcat('/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/IMERG_period/TEW_counter_clim2v2.mat');
% save(var_name3,'counter_clim2')

% Now ready to plot 2D column integrated avg LH for TEWs (1981-2018)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run in
%/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period
close all; clear;
load coast.mat
%load TEW_avgLHTEST2.mat
load TEW_sumLH.mat %the average of each yrly 2D integration
%load TEW_avgLH2.mat
LHavg_clim=LHavg_clim.*(1/9.81);
%avg_eqvRR=(LHavg_clim.*10.*1000./250).*365; %multiply by 365 to get ann avg accum in mm
avg_eqvRR=LHavg_clim.*10.*1000./250;



lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)'; 

c1=[1 1 1]; 
c2=[.94 1 1];
%c3=[.68 .85 .90];
c4=[.55 .75 .84];
%c5=[.27 .51 .71];
c6=[0 .65 .58];
%c7=[.19 .70 .10];
%c8=[.60 .80 .20];
c9=[.82 .89 .19];
c10=[1 .80 .20];
c11=[1 .47 0];
c12=[1 .22 0];
%c13=[.89 .26 .20];
%c14=[.77 .12 .23];
c15=[.62 .11 .20];
colors=[c1; c2; c4; c6; c9; c10; c11; c12; c15];
intervals=[0.1,1.0,5.0,10.0,50.0,100.0,250.0,500.0,1000.0];
normalized_intervals = (intervals - min(intervals)) / (max(intervals) - min(intervals));
colormap_custom = colormap(interp1(normalized_intervals,colors,linspace(0,1,9)));
color_indices=discretize(avg_eqvRR,intervals);

%GET GLOBAL MEAN CONDITIONAL MEAN
glob_mean = num2str(round((nansum(nansum(avg_eqvRR))))./...
    (sum(sum(~isnan(avg_eqvRR)))),1);

close all;
figure(1)
pcolor(lon_tdt,lat_tdt,(avg_eqvRR)'); shading interp;
colorbar;
colormap(gca,colormap_custom);
caxis([min(intervals) max(intervals)]);
hold on
plot(long,lat,'k','LineWidth',.75)
%shading flat
axis tight
axis equal
%cb=colorbar('horizonta','position',[0.2 0.3 0.6 0.01]); %left,bott,W,H
%cb.Label.String = 'mm'; %day^{-1}
axis([ -179.75 179.75 -36.75 36.75])
xticklabels({['150' char(176) 'W'],['100' char(176) 'W'],['50' char(176) 'W'],['0' char(176)],...
    ['50' char(176) 'E'],['100' char(176) 'E'],['150' char(176) 'E']})
yticklabels({['20' char(176) 'S'],['0' char(176)],['20' char(176) 'N']})
title(['                    Equivalent Surface Rain Rates (2001-2020)         \mu = ',glob_mean])
print('-f1','/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/accum_TEST1','-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear;
load coast.mat
% load TEW_avgLH_TEST.mat %should already be ann avg accum in mm
% %load TEW_avgLH2.mat
% LHavg_clim=LHavg_clim.*(1/9.81);
% avg_eqvRR=LHavg_clim.*10.*1000./250; 

load TEW_sumLH.mat %the average of each yrly 2D integration
%load TEW_avgLH2.mat
LHavg_clim=LH_accum_avg_clim./24.*3; %K/day to K/hr*3hrs
LHavg_clim=LHavg_clim.*(1/9.81);
%avg_eqvRR=(LHavg_clim.*10.*1000./250).*365; %multiply by 365 to get ann avg accum in mm
avg_eqvRR=LHavg_clim.*10.*1000./250;

%GET GLOBAL MEAN CONDITIONAL MEAN
glob_mean = num2str(round((nansum(nansum(avg_eqvRR))))./...
    (sum(sum(~isnan(avg_eqvRR)))),1);

c1=[1 1 1]; 
c2=[.94 1 1];
%c3=[.68 .85 .90];
c4=[.55 .75 .84];
%c5=[.27 .51 .71];
c6=[0 .65 .58];
%c7=[.19 .70 .10];
%c8=[.60 .80 .20];
c9=[.82 .89 .19];
c10=[1 .80 .20];
c11=[1 .47 0];
c12=[1 .22 0];
%c13=[.89 .26 .20];
%c14=[.77 .12 .23];
c15=[.62 .11 .20];
%mycolormap=rgbmap(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,15);
mycolormap=rgbmap(c1,c2,c4,c6,c9,c10,c11,c12,c15,9);



lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)'; 
close all;
figure(1)
pcolor(lon_tdt,lat_tdt,avg_eqvRR');
colormap(mycolormap)
caxis([0 1000]); set(gca,'ColorScale','log');
hold on
plot(long,lat,'k','LineWidth',.75)
shading flat
axis tight
axis equal
cb=colorbar('horizonta','position',[0.2 0.3 0.6 0.01]); %left,bott,W,H
cb.Label.String = 'mm'; %day^{-1}
axis([ -179.75 179.75 -36.75 36.75])
xticklabels({['150' char(176) 'W'],['100' char(176) 'W'],['50' char(176) 'W'],['0' char(176)],...
    ['50' char(176) 'E'],['100' char(176) 'E'],['150' char(176) 'E']})
yticklabels({['20' char(176) 'S'],['0' char(176)],['20' char(176) 'N']})
title(['                    Equivalent Surface Rain Rates (2001-2020)         \mu = ',glob_mean])
print('-f1','/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/annacumTEST_log','-dpng')




























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating domain means

%Conditional TEW LH rates

load TEW_avgLH.mat
LHavg_clim=LHavg_clim*(1/9.81); 

%GET GLOBAL MEAN
glob_mean = num2str(round((nansum(nansum(LHavg_clim)))./...
    (sum(sum(~isnan(LHavg_clim)))),2));
disp(strcat('Global mean LH rate (K/day) is _',glob_mean,' .'))

%GET NA MEAN 
%NA 5N-25N, 20W-70W; NEP 5N-25N, 90-140W
NA_box = LHavg_clim(177:257,91:131); %Draw box over N. Atlantic for TEW heating
NA_mean = num2str(round((nansum(nansum(NA_box)))./...
    (sum(sum(~isnan(NA_box)))),2));
disp(strcat('N.A. mean LH rate (K/day) is _',NA_mean,' .'))

%GET NEP MEAN 
NEP_box = LHavg_clim(65:145,91:131); %Draw box over N.E. Pacific for TEW heating
NEP_mean = num2str(round((nansum(nansum(NEP_box)))./...
    (sum(sum(~isnan(NEP_box)))),2));
disp(strcat('N.E.P. mean LH rate (K/day) is _',NEP_mean,' .'))



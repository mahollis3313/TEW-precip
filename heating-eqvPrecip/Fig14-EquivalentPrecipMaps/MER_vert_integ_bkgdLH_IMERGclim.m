
clear;

dirname1 = '/data/reanalysis/MERRA2/Daily/DTDT_TRMM_zlvls_v2/';
dirname2 = '/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/';
for y=2001:2020
    LH_vert_accum=zeros(576,161);
    counter_accum=zeros(576,161);
    for m=1:12
        days_in_mth = eomday(y,m); %get the number of days in each month
        for d=1:days_in_mth
            year = num2str(y); month = num2str(m); day = num2str(d);
            filename=sprintf('MERRA2_asm.tdt.3hr.%04d%02d%02d.nc', y, m, d); 
            fname = fullfile(dirname1, filename); %full filename with path
            LH_temp=ncread(fname,'LH_Vavg'); %already K/day
            LH_temp=squeeze(nansum(LH_temp,3));
            LH_vert_accum = LH_temp + LH_vert_accum;
            counter_temp = ~isnan(LH_temp); %assigns zeros to NaNs and 1's to non-NaNs
            counter_accum = counter_temp + counter_accum; 
        end
    end
    var_name = strcat(dirname2,'MERbkgdLHaccum',year,'.mat');
    save(var_name,'LH_vert_accum')
    var_name = strcat(dirname2,'MERbkgd_counter',year,'.mat');
    save(var_name,'counter_accum')
    LH_vert_rate = LH_vert_accum./counter_accum;
    var_name = strcat(dirname2,'MERbkgdLHrate',year,'.mat');
    save(var_name,'LH_vert_rate')
end

filelist1 = dir([dirname2 'MERbkgdLHaccum*.mat']);
filelist2 = dir([dirname2 'MERbkgd_counter*.mat']);
LH_accum = zeros(576,161,20); counter_accum = zeros(576,161,20);                     

for i = 1:20
    load(strcat(filelist1(i).name));
    LH_accum(:,:,i) = LH_vert_accum; %yearly accumulations of vert integ, LH (already includes 1/g)
    load(strcat(filelist2(i).name));
    counter_accum(:,:,i) = counter_accum;
end

bkgd_LHavg_accum = squeeze(nanmean(LH_accum,3)); %get annual average accumulation
bkgd_counter_accum = squeeze(nanmean(counter_accum,3)); %get annual average accumulated counts

var_name = strcat(dirname2,'bkgd_LHaccum_climavg.mat');
save(var_name,'bkgd_LHavg_accum')
var_name = strcat(dirname2,'bkgd_counter_climavg.mat');
save(var_name,'bkgd_counter_accum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate contributions from background (don't need to convert to mm to do this)

clear; close all;
%Get numerator (TEW or conditional LH accum) and denominator (unconditional bkgd accum)
%load TEW_avgLH.mat; TEW_LHavg_accum=LH_accum_avg_clim2; clear LH_accum_avg_clim; %the average of each yrly 2D integration
%load bkgd_LHaccum_climavg.mat; %the average of each yrly 2D integration
load TEW_avgLH.mat; load MERbkgdLHrateIMERG.mat;

%Load counter vars to get freq-weight
%load TEW_counter_climavg2.matl; load bkgd_counter_climavg.mat
load TEW_counter_clim.mat; load MERbkgd_counterIMERG.mat;

%freq_weight = TEW_counter_accum2 ./ bkgd_counter_accum; 
%TEW_contrib = (TEW_LHavg_accum./bkgd_LHavg_accum).*freq_weight; 
freq_weight = counter_clim ./ counter_accum;
TEW_contrib = (LHavg_clim./LH_vert_rate).*freq_weight;

%GET GLOBAL MEAN CONDITIONAL MEAN
glob_mean = num2str(round((nansum(nansum(TEW_contrib))))./...
    (sum(sum(~isnan(TEW_contrib)))),3);


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
pcolor(lon_tdt,lat_tdt,TEW_contrib'); 
hold on
plot(long,lat,'k','LineWidth',.75)
shading flat
axis tight
axis equal
colormap(mycolormap)
%rectangle('position',[-70 5 50.5 20.5 ],'EdgeColor','k','LineWidth',.75)
%rectangle('position',[-140 5 50.5 20.5 ],'EdgeColor','k','LineWidth',.75)
%rectangle('position',[-140 -25 50.5 20.5 ],'EdgeColor','k','LineWidth',.75)
%rectangle('position',[-50 -25 61 20.5 ],'EdgeColor','k','LineWidth',.75)
caxis([0 .15]);
cb=colorbar('horizonta','position',[0.2 0.3 0.6 0.01]); %left,bott,W,H
%cb.Label.String = 'K day^{-1}';
axis([ -179.75 179.75 -36.75 36.75])
xticklabels({['150' char(176) 'W'],['100' char(176) 'W'],['50' char(176) 'W'],['0' char(176)],...
    ['50' char(176) 'E'],['100' char(176) 'E'],['150' char(176) 'E']})
yticklabels({['20' char(176) 'S'],['0' char(176)],['20' char(176) 'N']})
title(['                    Fractional Contribution from TEWs (2001-2020)     \mu = ',glob_mean])
print('-f1','MER_Vinteg_fractional_contrib_IMERG_850waves-cb4-weightedv2-updated','-dpng')




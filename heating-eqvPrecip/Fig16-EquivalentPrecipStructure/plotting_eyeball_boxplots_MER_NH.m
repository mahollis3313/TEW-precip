%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% MER NH TEW Accum/Avg LH (2001-2020) %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear; close all; 

folder = '/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/eyeballs700/';
filelist1 = dir([folder 'NH_TEW_accumavgLH_*.mat']);
filelist2 = dir([folder 'NH_TEW_latlon_counter_*.mat']);

LH_ann = zeros(17,21,20);                             
counter_ann = zeros(17,21,20); 
for i = 1:20
    load(strcat(folder,(filelist1(i).name)));
    LH_ann(:,:,i) = TEW_integ_heating; %K/day
    load(strcat(folder,(filelist2(i).name)));
    counter_ann(:,:,i) = counter_integ; %counter for all waves
end

LH_accum_avg_clim = squeeze(nanmean(LH_ann,3));
counter_clim = squeeze(nanmean(counter_ann,3)); %should this use sum instead of mean? TEST WITH BOTH
LHavg_clim = LH_accum_avg_clim ./ counter_clim;

counter_clim_all = squeeze(nanmean(counter_ann,3)); 

var_name = strcat(folder,'NH_TEW_accumavgLH.mat');
save(var_name,'LH_accum_avg_clim')

var_name1 = strcat(folder,'NH_TEW_avgLH.mat');
save(var_name1,'LHavg_clim')

var_name2 = strcat(folder,'NH_TEW_counter_ALL.mat');
save(var_name2,'counter_clim_all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%For all years (1981-2018) - plotting conditional AVERAGE TEW LH for N.H.
close all; clear;
folder = '/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/eyeballs700/';
load(strcat(folder,'NH_TEW_avgLH.mat'));
LHavg_clim=LHavg_clim.*(1/9.81); %add in 1/g multiplier

LHavg_clim(LHavg_clim==0)=NaN; 
lat_MER=1:21; lon_MER=1:17;

close all;
figure(1)
%pcolor(lon_MER,lat_MER,LHavg_clim); 
[C,h]=contourf((LHavg_clim.*40./24)',15);  %*40 to take to mm/day then /24 to take to mm/hr
w=h.LevelStep;
h.LevelStep = 0.5; 
hold on
shading flat
axis tight
axis square
Ax = gca; Ax.YGrid = 'on'; Ax.XGrid = 'on'; Ax.GridAlpha = 1.0; Ax.Layer = 'top'; Ax.GridColor = [0,0,0];
pbaspect([10 10 1]); hold on;
% colormap(parula(12))
% cbh.Ticks = [0 1 2 3 4 5 6 7 8 12];
% caxis([0 12]);
colormap(parula(15))
%cbh.Ticks = [0 1 2 3 4 5 6 7 8 12];
caxis([0 1.5]);
colorbar
h = colorbar; ylabel(h,'mm hr^{-1}')%(h,'Q_{1} (K day^{-1})');
xlim([1 17])
xticks([1 3 5 7 9 11 13 15 17])
xticklabels({'-5','-3.125','-2.5','-1.25','0','1.25','2.5','-3.125','5'})
%ylim([3 19])
yticks([1 3 5 7 9 11 13 15 17 19 21])
yticklabels({'-5','-4','-3','-2','-1','0','1','2','3','4','5'})
xlabel('Degrees from center')
ylabel('Degrees from center')
title({'Equivalent Surface Rain Rates (2001-2020)', 'within 500-km of TEW kinematic center'})
print('-f1','NH_700waves_eqvrain_eyeball_v2-15int','-dpng')



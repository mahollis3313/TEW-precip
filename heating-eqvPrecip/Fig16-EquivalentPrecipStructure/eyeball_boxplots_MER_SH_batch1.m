%BATCH 1 of 2

%3-hourly 850 hPa wave tracks for S. HEMISPHERE
%matching to MERRA-2 'LH_Vavg' to get accumulated rates in N.H.
%VERTICALLY INTEGRATED - RESOLVED FOR 3hrly matches and not daily integral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear;
lat_tdt=(-40:.5:40)'; lon_tdt=(-180:.625:179.3750)'; wave_count = 0;

dirname = '/home/c770m986/TEW_MER_tracking/FilteredNoTCUpdated/700hPa/'; %change wave database here 850 or 700
for y = 2001:2020
    year = num2str(y);
    Strack_fname=strcat(dirname,'post_filter.',year,'.SH.NoTC.nc');
    
    %Read in the tracking info
    lon=ncread(Strack_fname,'longitude'); lon=wrapTo180(lon);
    lat=ncread(Strack_fname,'latitude');
    time=ncread(Strack_fname,'time'); wave_count = wave_count + length(time);

    time=datetime(1979,12,1,0,0,0) + hours(time);
    time.Format='yyyyMMddHHmmss';  
    formatOut='yyyymmddHHMMss'; 
    time=datestr(time,formatOut);

    %Initialize 
    TEW_integ_heating=zeros(17,21); %x,y or (17 x 21 x 36)
    counter_integ = zeros(17,21);
   
    %First temporal matching
    for i = 1:length(time) 
        %if (lat(i) < 25.25) && (lat(i) > 4.75) 
            %if (lon(i) > -70.3125) && (lon(i) < -19.6875)  %+/- 1/2*resolution
                t_w = time(i,:);
                f_w = t_w(1:8);
                f_list = strcat('/data/reanalysis/MERRA2/Daily/DTDT_TRMM_zlvls_v2/MERRA2_asm.tdt.3hr.',f_w,'.nc');
                dtdt_time = ncread(f_list,'time');
                yr = t_w(9:10);
                yr_out = str2double(yr);
                t_mm = yr_out*60; %convert from UTC hour to UTC minutes
                t_index = find(dtdt_time==t_mm);
                LH_Vavg = ncread(f_list,'LH_Vavg');  
                
                %spatial matching at center -/+ 5 degrees 
                %lon matching using physical degrees not indices
%                 lon_window1=lon(i)-5; lon_window2=lon(i)+5; 
%                 lon_diff1 = abs(lon_tdt - lon_window1); 
%                 lon_min1 = min(lon_diff1);
%                 lon_index1 = find(lon_diff1==lon_min1);
%                 lon_diff2 = abs(lon_tdt - lon_window2); 
%                 lon_min2 = min(lon_diff2);
%                 lon_index2 = find(lon_diff2==lon_min2);
                
                lon_diff = abs(lon_tdt - lon(i)); 
                lon_min = min(lon_diff);
                lon_index = find(lon_diff==lon_min);
                if size(lon_index,1) > 1
                    lon_index=lon_index(1); 
                else
                    %do nothing
                end
                lon_window1=lon_index-8; lon_window2=lon_index+8; %8 indices = 5 degrees lon
                
                if lon_window2 > 576 %if lon index > 576
                    temp_lon1 = 1; temp_lon2 = 17-(lon_window2-576);
                    lon_range = lon_window1:576;
                elseif lon_window1 < 1 %if lon index < 1
                    temp_lon1 = -(lon_window1-1); temp_lon2 = 17;
                    lon_range = 1:lon_window2+1;   %should it be +1 here?
                elseif lon_window1 <= 576 && lon_window2 >= 1
                    temp_lon1 = 1; temp_lon2 = 17;
                    lon_range = lon_window1:lon_window2;
                end
                
                %lat matching using indices not physical degrees
                lat_diff = abs(lat_tdt - lat(i)); 
                lat_min = min(lat_diff);
                lat_index = find(lat_diff==lat_min);
                if size(lat_index,1) > 1
                    lat_index=lat_index(1); 
                else
                    %do nothing
                end
                lat_window1=lat_index-10; lat_window2=lat_index+10; %10 indices = 5 degrees lat
                
                if lat_window2 > 161 %&& lat_window1 >= 1 %if lat is outside of 40 N
                    temp_lat1 = 1; temp_lat2 = 21-(lat_window2-161);
                    lat_range = lat_window1:161;
                elseif lat_window1 < 1 %&& lat_window2 <= 576 %if lat is outside of 40 S
                    temp_lat1 = -(lat_window1-1); temp_lat2 = 21;
                    lat_range = 1:lat_window2+1;
                elseif lat_window1 <= 161 && lat_window2 >= 1
                    temp_lat1 = 1; temp_lat2 = 21;
                    lat_range = lat_window1:lat_window2;
                end        

                
                temp_heating = nan(17,21);
                temp_heating(temp_lon1:temp_lon2,temp_lat1:temp_lat2) =...
                    LH_Vavg(lon_range(1):lon_range(length(lon_range)),lat_range(1):lat_range(length(lat_range)),t_index);
               % temp_heating(:,temp_lat1:temp_lat2) = LH_Vavg(lon_index1:lon_index2,lat_range(1):lat_range(length(lat_range)));
                
                temp_counter = zeros(17,21);
                temp_counter(~isnan(temp_heating)) = 1;
                counter_integ = counter_integ + temp_counter; %'accumulated' counts (17x21)
                temp_heating(isnan(temp_heating)) = 0; %replace all NaN values with zeros, while keeping all other values the same
                TEW_integ_heating = TEW_integ_heating + temp_heating; %accumlated heating (17x21)
                
                %clear temp_lon; clear lon_range; %shouldn't need clear statements?
                %clear temp_lat; clear lat_range;
                
            %else
                %do nothing since lon falls outside of NA
            %end
        %else
            %do nothing since lat falls outside of domain
        %end
        
      frac_now =5*(ceil(i/length(time)*100/5));
      frac_next=5*(ceil((i+1)/length(time)*100/5));
      if frac_next>frac_now; disp([num2str(frac_now) '% completed']); end 
    
    end
    
    
    %already in K/day?
    %TEW_integ_heating = TEW_integ_heating*86400; %Leave as accumulated rate (K/day)          *(1/8); %K/sec to K/day to K


folder='/home/c770m986/MATLAB_scripts/a_manuscript_code/new_analysis/MER/fixingMER_VintegLH/IMERG_period/eyeballs700/';
var_name1 = strcat(folder,'SH_TEW_accumavgLH_700_',year,'.mat');
save(var_name1,'TEW_integ_heating')               

var_name2 = strcat(folder,'SH_TEW_latlon_counter_700_',year,'.mat');
save(var_name2,'counter_integ')     

end

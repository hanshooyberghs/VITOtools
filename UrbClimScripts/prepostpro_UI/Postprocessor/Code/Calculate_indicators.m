function results_indicators=Calculate_indicators(daily_data,indicator_data)

list_indicators=fieldnames(indicator_data);

for i=1:length(list_indicators)
    indicator=list_indicators{i};
    if indicator_data.(indicator).use ==1
        disp(indicator)
        
        switch indicator
            case 'Tmin_mean'
                results_indicators.Tmin_mean.data=mean(daily_data.Tmin,3);
                results_indicators.Tmin_mean.name='Tmin_mean';
            case 'Tmean_mean'
                results_indicators.Tmean_mean.data=mean(daily_data.Tmean,3);
                results_indicators.Tmean_mean.name='Tmean';
            case 'Tmax_mean'
                results_indicators.Tmax_mean.data=mean(daily_data.Tmax,3);
                results_indicators.Tmax_mean.name='Tmax_mean';
                
            case 'Tmin_prctile'
                results_indicators.Tmin_prctile.data=prctile(daily_data.Tmin,indicator_data.Tmin_prctile.prctile,3);
                results_indicators.Tmin_prctile.name=['Tmin_p',num2str(indicator_data.Tmin_prctile.prctile)];
            case 'Tmean_prctile'
                results_indicators.Tmean_prctile.data=prctile(daily_data.Tmean,indicator_data.Tmean_prctile.prctile,3);
                results_indicators.Tmean_prctile.name=['Tmean_p',num2str(indicator_data.Tmean_prctile.prctile)];
            case 'Tmax_prctile'
                results_indicators.Tmax_prctile.data=prctile(daily_data.Tmax,indicator_data.Tmax_prctile.prctile,3);
                results_indicators.Tmax_prctile.name=['Tmax_p',num2str(indicator_data.Tmax_prctile.prctile)];
                
            case 'Tfixedhour_mean'
                hour=indicator_data.Tfixedhour_mean.hour;
                results_indicators.Tfixedhour_mean.data=mean(daily_data.selectedhour,3);
                results_indicators.Tfixedhour_mean.name=['T_',num2str(hour),'hUTC'];
                
            case 'cdd_cumul'
                threshold=indicator_data.cdd_cumul.temperature;
                results_indicators.cdd_cumul.data=sum(max(daily_data.Tmax-threshold,0),3);
                results_indicators.cdd_cumul.name=['CoolingDegreeDays_',num2str(threshold),'degrees'];
                
            case 'cdd_cumul_mean'
                threshold=indicator_data.cdd_cumul_mean.temperature;
                results_indicators.cdd_cumul_mean.data=sum(max(daily_data.Tmean-threshold,0),3);
                results_indicators.cdd_cumul_mean.name=['CoolingDegreeDaysMean_',num2str(threshold),'degrees'];
            
            case 'trop_night'
                threshold=indicator_data.trop_night.temperature;
                results_indicators.trop_night.data=sum(daily_data.Tmin>threshold,3);
                results_indicators.trop_night.name=['TropicalNights_above',num2str(threshold),'degrees'];
                
            case 'trop_day'
                threshold=indicator_data.trop_day.temperature;
                results_indicators.trop_day.data=sum(daily_data.Tmax>threshold,3);
                results_indicators.trop_day.name=['TropicalDays_above',num2str(threshold),'degrees'];
                
            case 'trop_day_mean'
                threshold=indicator_data.trop_day_mean.temperature;
                results_indicators.trop_day_mean.data=sum(daily_data.Tmean>threshold,3);
                results_indicators.trop_day_mean.name=['TropicalDaysMean_above',num2str(threshold),'degrees'];
                
            case 'hwd'
                results_indicators.hwd=calc_heatwaves(daily_data,indicator_data.hwd);
                                           
        end                   
    end    
end


end

function hwd = calc_heatwaves(daily_data,hwd_input)
    if hwd_input.runningmean==1
    %% option 1: running mean
        %calculate running means
        Tmin_running=movingmean(daily_data.Tmin,hwd_input.nrdays,3);
        Tmax_running=movingmean(daily_data.Tmax,hwd_input.nrdays,3);

        temp=(Tmin_running>hwd_input.tmin)+(Tmax_running>hwd_input.tmax);
        hwdres=sum(temp>1.5,3);
        name=['HeatWaveDays_RunningThresholds_Tmin',num2str(hwd_input.tmin),'Tmax_',num2str(hwd_input.tmax)];

    else
    %% option 2: no running mean (look at periods with exceedances)
        hwdres=zeros(size(daily_data.Tmin));

        % find last day of period    
        locdag2=intersect(find(daily_data.Tmin>hwd_input.tmin),find(daily_data.Tmax>hwd_input.tmax));

        % check previous days for all locations
        for i=1:length(locdag2)
           [x,y,day]=ind2sub(size(daily_data.Tmin),locdag2(i));       
           % only look at days that are not at the start of the period
           if day - hwd_input.nrdays>0
               yes=1;
               j=1;
               while j<hwd_input.nrdays
                   if daily_data.Tmax(x,y,day-j)<hwd_input.tmax
                       yes =0;
                       j=100;
                   elseif daily_data.Tmin(x,y,day-j)<hwd_input.tmin
                       yes =0;
                       j=100;
                   end
                   j=j+1;
               end

               if yes ==1
                    hwdres(x,y,day-hwd_input.nrdays+1:day)=1;
               end
           end
        end
        hwdres=sum(hwdres,3);
        name=['HeatWaveDays_Thresholds_Tmin',num2str(hwd_input.tmin),'Tmax_',num2str(hwd_input.tmax)];

    end
    hwd.data=hwdres;
    hwd.name=name;
        
end
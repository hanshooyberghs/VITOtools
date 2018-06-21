function indicator_data= SetUpIndicators()

%create default indicator values.
clear('indicator_data')

% Tmin mean
indicator_data.Tmin_mean.use=1;

% Tmax mean
indicator_data.Tmean_mean.use=0;

% Tmax_mean
indicator_data.Tmax_mean.use=0;

% Tmin mean
indicator_data.Tmin_prctile.use=0;
indicator_data.Tmin_prctile.prctile=95;

% Tmax mean
indicator_data.Tmean_prctile.use=0;
indicator_data.Tmean_prctile.prctile=95;

% Tmax_mean
indicator_data.Tmax_prctile.use=0;
indicator_data.Tmax_prctile.prctile=95;

% Temperature at XXh UTC
indicator_data.Tfixedhour_mean.use=1;
indicator_data.Tfixedhour_mean.hour=2;

% Cooling degree cumulative 
indicator_data.cdd_cumul.use=0;
indicator_data.cdd_cumul.temperature=25;

% Tropical nights 
indicator_data.trop_night.use=0;
indicator_data.trop_night.temperature=18;

% Tropical nights 
indicator_data.trop_day.use=0;
indicator_data.trop_day.temperature=30;


% Cooling degree cumulative  mean
indicator_data.cdd_cumul_mean.use=0;
indicator_data.cdd_cumul_mean.temperature=25;


% Tropical nights 
indicator_data.trop_day_mean.use=0;
indicator_data.trop_day_mean.temperature=30;


% HWD
indicator_data.hwd.use=0;
indicator_data.hwd.nrdays=3;
indicator_data.hwd.tmin=18;
indicator_data.hwd.tmax=30;
indicator_data.hwd.runningmean=1;

end
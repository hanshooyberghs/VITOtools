%% main routines
function [results_indicators,data_geo]=ProcessUrbClim(file_in,indicator_data,output_formats)

disp('Read data')
if iscell(file_in)
% case with multiple input files
    nr_files=length(file_in)-1;
    for i=2:length(file_in)
        file=[file_in{1},file_in{i}];
        disp(['Read file: ',file])
        [daily_data_file,data_geo_file]=verwerk_file(file,indicator_data.Tfixedhour_mean.hour);
        
        % process data
        if i==2 % first run
            data_geo=data_geo_file;
            daily_data=daily_data_file;
        else
            % check if geo data matches previous data
            if ~isequal(data_geo,data_geo_file)
                err = MException('InputError:InputGeo', ...
                    'Spatial data of files does not match');
                throw(err)
            end
            
            % add data to daily data
            daily_data.date=cat(1,daily_data.date,daily_data_file.date);
            daily_data.Tmin=cat(3,daily_data.Tmin,daily_data_file.Tmin);
            daily_data.Tmean=cat(3,daily_data.Tmean,daily_data_file.Tmean);
            daily_data.Tmax=cat(3,daily_data.Tmax,daily_data_file.Tmax);
            daily_data.selectedhour=cat(3,daily_data.selectedhour,daily_data_file.selectedhour);
        end
        
    end

else
% case with only one input file
    disp(['Read one file: ',file_in])
    nr_files=1;
    [daily_data,data_geo]=verwerk_file(file_in,indicator_data.Tfixedhour_mean.hour);
end

fprintf('\n\nCalculate indicators\n')
results_indicators=Calculate_indicators(daily_data,indicator_data);
if isfield(results_indicators, 'hwd')
    disp(nr_files)
    results_indicators.hwd.data=results_indicators.hwd.data/nr_files;
end
fprintf('\n\nWrite to Output\n')
write_data(results_indicators,data_geo,output_formats)


end

%% Process one file to geo-data and daily temperature data
function [daily_data,data_geo]=verwerk_file(file_in,hour_mean)
    [data_temperature,data_geo]=read_file(file_in);
    daily_data=convert_to_daily(data_temperature);
    daily_data.selectedhour=data_temperature.t2m(:,:,data_temperature.hour==hour_mean);
    clear data_temperature
end



%% reading input data
function [data_temperature,data_geo]=read_file(file_in)
%% read temperature
t_in=ncread(file_in,'Surface/T2M')-273.15;
temp=size(t_in);
grootte=temp(3);
t_in=t_in(:,:,1:grootte-1);  %to remove for use with new urbclim

%% read time
time_in=ncread(file_in,'Parameters/Time');
time_in=time_in(:,1:grootte-1);   %to remove for use with new urbclim
year=str2num(time_in(1:4,:)');
month=str2num(time_in(6:7,:)');
day=str2num(time_in(9:10,:)');
hour=str2num(time_in(12:13,:)');
date=datenum(year,month,day,hour,zeros(size(hour)),zeros(size(hour)));


%% read geographic data
x_in=ncread(file_in,'Parameters/X');
y_in=ncread(file_in,'Parameters/Y');
lat_in=ncread(file_in,'Parameters/LAT');
lon_in=ncread(file_in,'Parameters/LON');



ncid = netcdf.open(file_in,'NC_NOWRITE');
espg_in=str2num(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'EPSG'));
netcdf.close(ncid);

%% group in structure
data_temperature.t2m=t_in; 
data_temperature.year=year;
data_temperature.month=month;
data_temperature.day=day;
data_temperature.hour=hour;
data_temperature.date=date;
data_geo.x=x_in;
data_geo.y=y_in;
data_geo.lat=lat_in;
data_geo.lon=lon_in;
data_geo.EPSG=espg_in;

end

%% find daiily min, mean and max values
function daily_data=convert_to_daily(data_temperature)
temp=size(data_temperature.t2m);
nx=temp(1);
ny=temp(2);

list_days=unique(floor(data_temperature.date));
maxdays=zeros(nx,ny,length(list_days));
mindays=zeros(nx,ny,length(list_days));
gemdays=zeros(nx,ny,length(list_days));

for i=1:length(list_days)
    day=list_days(i);
    maxdays(:,:,i)=max(data_temperature.t2m(:,:,floor(data_temperature.date)==day),[],3);
    mindays(:,:,i)=min(data_temperature.t2m(:,:,floor(data_temperature.date)==day),[],3);
    gemdays(:,:,i)=mean(data_temperature.t2m(:,:,floor(data_temperature.date)==day),3);
end
daily_data.date=list_days;
daily_data.Tmin=mindays;
daily_data.Tmax=maxdays;
daily_data.Tmean=gemdays;


end


%% write outputdata
function write_data(results_indicators,data_geo,output_formats)
    mkdir('Results')
    
    % set colormap data
    load('colordata.mat')
    figure('visible','off')
    colormap(colordata)
    close(gcf)


    resolutie_x=max(diff(data_geo.x(:)));
    resolutie_y=max(diff(data_geo.y(:)));
    epsg=data_geo.EPSG;

    list_indicators=fieldnames(results_indicators);
    for i=1:length(list_indicators)
        % get indicator results
        indicator=list_indicators{i};
        disp(indicator)
        data=results_indicators.(indicator).data;
        name=results_indicators.(indicator).name;
        
        %save matlab
        save(['Results/',name,'.mat'],'data');
        
        
        % deduce plotting vector
        min_dat=min(data(:));
        max_dat=max(data(:));
        step=(max_dat-min_dat)/30;
        vector=[min_dat:step:max_dat];        
        
        % save image
        if output_formats.png==1
            h=m_map_pcolor(data_geo.lon,data_geo.lat,data,colordata);
            outputfile=['Results/',name,'.png'];
            print(h,'-dpng','-r300',outputfile)
            set(h,'PaperUnits','inches','PaperPosition',[0 0 30 30])
            set(findall(h,'-property','FontSize'),'FontSize',25)
            close(h)       
        end
        
        % make kmz 
        if output_formats.kmz==1
            writeKMZ(data_geo.lon,data_geo.lat,data,colordata,name)
        end
        
        % define geotiff structure and write to output file
        if output_formats.tif==1
            outputfile=['Results/',name,'.tif'];
            geo = struct('ModelPixelScaleTag', [ resolutie_x resolutie_y 0 ], ...
                'GTModelTypeGeoKey',1, ...
                'ProjectedCSTypeGeoKey',epsg ,... 
                'ModelTiepointTag',   [ 0.5 0.5 0 min(data_geo.x(:)) max(data_geo.y(:)) 0 ]);   
            geotiffwrite( outputfile, [], flipud(data'), 32, geo );  
        end
    
    end
end

function writeKMZ(lon,lat,data,colordata,name)

    temp=size(lat);
    nx=temp(1);
    ny=temp(2);
    a=kml(name);

    increment=(max(data(:))-min(data(:)))/(length(colordata)-1);
    splits=[min(data(:)):increment: max(data(:))];
    splits=[splits max(data(:))+1];

    for i=2:nx-1
        for j=2:ny-1
            temp=find(data(i,j)-splits<0);
            colorHex= kml.color2kmlHex([ colordata(temp(1)-1,:) 1] );
            lon1=(lon(i,j)+lon(i-1,j-1))/2;
            lat1=(lat(i,j)+lat(i-1,j-1))/2;
            lon2=(lon(i,j)+lon(i+1,j-1))/2;
            lat2=(lat(i,j)+lat(i+1,j-1))/2;
            lon3=(lon(i,j)+lon(i+1,j+1))/2;
            lat3=(lat(i,j)+lat(i+1,j+1))/2;
            lon4=(lon(i,j)+lon(i-1,j+1))/2;
            lat4=(lat(i,j)+lat(i-1,j+1))/2;
            subarray_lon=[lon1 lon2 lon3 lon4];
            subarray_lat=[lat1 lat2 lat3 lat4];      
            a.poly3(subarray_lon,subarray_lat,zeros(size(subarray_lat)), 'polyColor', colorHex, ...
               'altitudeMode','clampToGround', 'lineWidth',0,'name',['poly_' num2str(i) num2str(j)], ...
               'id',[ 'poly_' num2str(i) num2str(j)] );

        end
    end
    disp('Start saving KMZ')
    outputfile=['Results/',name,'.kmz'];
    a.save(outputfile)
end



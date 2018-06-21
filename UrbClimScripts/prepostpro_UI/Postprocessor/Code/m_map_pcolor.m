function h=m_map_contourf(lon,lat,data,colordata)

%% function to plot data on a geographic map

% geographic processing
lon=double(lon);
lat=double(lat);
data=double(data);

min_lon=min(lon(:));
max_lon=max(lon(:));

min_lat=min(lat(:));
max_lat=max(lat(:));



%% plotting
m_proj('lambert','long',[min_lon max_lon],'lat',[min_lat max_lat]);

if min(size(lon))==1 && min(size(lat))==1,
 [lon,lat]=meshgrid(lon,lat);
end;
 
[X,Y]=m_ll2xy(lon,lat,'clip','on');  %First find the points outside
 
i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN; 
if any(i(:)), [X,Y]=m_ll2xy(lon,lat,'clip','patch'); end;  

if any(~i(:)),
   h=figure('visible','off'); 
   p=pcolor(X,Y,data); 
   set(p,'EdgeColor', 'none');
   colorbar
else
  cs=[];h=[];
end;


%% nice lay-out
m_grid('box','fancy','tickdir','in');
colormap(colordata)

end
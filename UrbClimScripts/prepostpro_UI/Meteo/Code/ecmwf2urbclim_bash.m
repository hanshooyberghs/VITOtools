function ecmwf2urbclim_bash(config_file)
%
% Authors original IDL-version: Dirk Lauwaet (dirk.lauwaetATvito.be) and Koen De Ridder
% (koen.deridderATvito.be) 
% Author MatlabVersion: Hans Hooyberghs (hans.hooyberghsATvito.be)
% VITO (c) 2014
%
% PROGRAM SUMMARY
% ---------------
% This code interpolates vertical profiles of relevant meteorological quantities 
% from ECMWF output fields to the center of the URBCLIM domain
% Input variables:
%       config files:   Traditional UrbClim config file
%
%
% Input data used:
% * ECMWF 3D fields
%          u,v:  horizontal wind speed components
%          t:    air temperature
%          q:    specific humidity
% * ECMWF 2D fields
%          sp:   surface pressure 
%          sshf: surface sensible heat flux
%          slhf: surface latent heat flux
%          ewss: east-west surface stress
%          nsss: north-south surface stress
%          ssrd: downwelling short-wave radiation
%          strd: downwelling long-wave radiation
%          swvlx:volumetric soil water in layer x
%          stlx: soil temperature in layer x
%          cp:   convective precipitation
%          lsp:  large-scale precipitation
%          z:    surface height
%          lsm:  land-sea mask
%          sst:  sea surface temperature
%        
% REMARK: in the data arrays
% - longitude goes from west to east
% - latitude goes from north to south
% - vertical levels go from high to low
%
%% process input => read config file
fid=fopen(config_file);
data=textscan(fid,'%s ','Delimiter','\n');
for i=1:length(data{1})
   temp=strsplit(data{1}{i},'=');
   switch strtrim(temp{1})
       case 'fname_ecmwf'
           outputfile=strtrim(temp{2});
           disp(['fname_ecmwf: ',outputfile])
           
      case 'fname_atmo2D'
           inputfile2d=strtrim(temp{2});
           disp(['fname_atmo2D: ',inputfile2d])
           
     case 'fname_atmo3D'
           inputfile3d=strtrim(temp{2});
           disp(['fname_atmo3D: ',inputfile3d])

     case 'fname_parms'
           modellevels_file=strtrim(temp{2});
           disp(['fname_parms: ',modellevels_file])
           
      case 'nz'
           nr_vertlevels=str2double(temp{2});
           disp(['nz: ',num2str(clat)])
           
       case 'clon'
           clon=str2double(temp{2});
           disp(['clon: ',num2str(clon)])
           
       case 'clat'
           clat=str2double(temp{2});
           disp(['clat: ',num2str(clat)])
           
       case 'tmet'
           tmet=str2double(temp{2});    
           disp(['tmet: ',num2str(tmet)])
   end
end

%% fixed input (to change in case you want special things)
%  logfolder
    splitted_name=strsplit(outputfile,'/');
    pretext=strrep(splitted_name{length(splitted_name)},'.nc','');
    logfolder='log/';
%  soil texture parameters (see config file)
    etas = 0.395;
    psis = 0.121;
    bch = 4.05;
    
% reference date
    refdate=datenum(strrep(ncreadatt(inputfile2d,'/time','units'),'hours since ',''))
    
% physical constants 
    grav = 9.80665;
    Rdry = 287.04;
    Rvap = 461.50;
    re   = 6371000;
    rot  = 7.2921e-5; %Earth's rotation rate
    mkdir(logfolder)
    
% adapt vertical levels    
    zout=[10.,35.,70.,120.,175.,240.,320.,410.,520.,650.,800.,1000.,1250.,1500.,1750.,2000.,2250.,2500.,2750.,3000];
    zout=zout(1:nr_vertlevels);

    

%% read data
% inlezen tijd
time=ncread(inputfile2d,'time');
timebis=ncread(inputfile3d,'time');

% check time formats
if ~isequal(time,timebis)    
    err = MException('InputError:InputECMWF', ...
        'Time format of 2d and 3d file do not match');
    throw(err)
end
level=ncread(inputfile3d,'level');
nz=length(level);   % aantallevels

time=double(time);
matlabtime=refdate+time/24;
ntime=length(matlabtime);

% inlezen 2d ecmwf input
lon=ncread(inputfile2d,'longitude');
lat=ncread(inputfile2d,'latitude');
ps=ncread(inputfile2d,'sp');
hh=ncread(inputfile2d,'sshf');
ee=ncread(inputfile2d,'slhf');
taux=ncread(inputfile2d,'ewss');
tauy=ncread(inputfile2d,'nsss');
ssrd=ncread(inputfile2d,'ssrd');
strd=ncread(inputfile2d,'strd');
wsl1=ncread(inputfile2d,'swvl1');
wsl2=ncread(inputfile2d,'swvl2');
wsl3=ncread(inputfile2d,'swvl3');
wsl4=ncread(inputfile2d,'swvl4');
tsl1=ncread(inputfile2d,'stl1');
tsl2=ncread(inputfile2d,'stl2');
tsl3=ncread(inputfile2d,'stl3');
tsl4=ncread(inputfile2d,'stl4');
cp=ncread(inputfile2d,'cp');
lsp=ncread(inputfile2d,'lsp');
gs=ncread(inputfile2d,'z');
esst=ncread(inputfile2d,'sst');
lsm=ncread(inputfile2d,'lsm');


% inlezen 3d ecmwf input
levels=ncread(inputfile3d,'level');
nlevs=length(levels);
qq=ncread(inputfile3d,'q');
tt=ncread(inputfile3d,'t');
uu=ncread(inputfile3d,'u');
vv=ncread(inputfile3d,'v');

if mean(lon) >180
    lon = lon-360;
end



%% check on longitude
if (clon > max(lon(:))) || (clon <min(lon(:)))
        err = MException('InputError:CoordinatesError', ...
        'Coordinates do not match with input file.');
    throw(err)
end

if (clat > max(lat(:))) || (clat <min(lat(:)))
        err = MException('InputError:CoordinatesError', ...
        'Coordinates do not match with input file.');
    throw(err)
end

%% basic verwerking
% corrigeer error in z
for i=1:ntime
 gs(:,:,i)=gs(:,:,1);
end
zsfc=gs(:,:,1)/9.8065;  %geopotential -> meters


%% recalculations
% some conversions
hh=-hh/10800;
ee=-ee/10800;
ssrd=ssrd/10800;
strd=strd/10800;
taux=taux/10800;
tauy=tauy/10800;
lsp=lsp*1000;
cp=cp*1000;

%convert to instantaneous valous
cpmod=zeros(size(cp));
lspmod=zeros(size(lsp));
ssrdmod=zeros(size(ssrd));
strdmod=zeros(size(strd));
eemod=zeros(size(ee));
hhmod=zeros(size(hh));
tauxmod=zeros(size(taux));
tauymod=zeros(size(tauy));

%% convert to instantaneous values
% hours before 3 or 15 => use first step at 3 or 15 (based on previous
% hours, so no correct value possible)

% find first occurence of 3 or 15 UTC

datevector=datevec(matlabtime);
hour=datevector(:,4);
loc_3_15=find(mod(hour,12)==3);
first_time=loc_3_15(1);

% replace everything up to that hour with the value for 3 o clock
cpmod(:,:,1:first_time-1)=repmat(cp(:,:,first_time),1,1,first_time-1);
lspmod(:,:,1:first_time-1)=repmat(lsp(:,:,first_time),1,1,first_time-1);
ssrdmod(:,:,1:first_time-1)=repmat(ssrd(:,:,first_time),1,1,first_time-1);
strdmod(:,:,1:first_time-1)=repmat(strd(:,:,first_time),1,1,first_time-1);
eemod(:,:,1:first_time-1)=repmat(ee(:,:,first_time),1,1,first_time-1);
hhmod(:,:,1:first_time-1)=repmat(hh(:,:,first_time),1,1,first_time-1);
tauxmod(:,:,1:first_time-1)=repmat(taux(:,:,first_time),1,1,first_time-1);
tauymod(:,:,1:first_time-1)=repmat(tauy(:,:,first_time),1,1,first_time-1);    


% for the moment: assume start at 3 or 15 UTC
for cnt=first_time: length(matlabtime)
   if mod(cnt-first_time,4)==0
       cpmod(:,:,cnt)=cp(:,:,cnt);
       lspmod(:,:,cnt)=lsp(:,:,cnt);
       ssrdmod(:,:,cnt)=ssrd(:,:,cnt);
       strdmod(:,:,cnt)=strd(:,:,cnt);
       eemod(:,:,cnt)=ee(:,:,cnt);
       hhmod(:,:,cnt)=hh(:,:,cnt);
       tauxmod(:,:,cnt)=taux(:,:,cnt);
       tauymod(:,:,cnt)=tauy(:,:,cnt);    
   else
       cpmod(:,:,cnt)=cp(:,:,cnt)-cp(:,:,cnt-1);
       lspmod(:,:,cnt)=lsp(:,:,cnt)-lsp(:,:,cnt-1);
       ssrdmod(:,:,cnt)=ssrd(:,:,cnt)-ssrd(:,:,cnt-1); 
       strdmod(:,:,cnt)=strd(:,:,cnt)-strd(:,:,cnt-1);    
       hhmod(:,:,cnt)=hh(:,:,cnt)-hh(:,:,cnt-1);    
       eemod(:,:,cnt)=ee(:,:,cnt)-ee(:,:,cnt-1);    
       tauxmod(:,:,cnt)=taux(:,:,cnt)-taux(:,:,cnt-1);
       tauymod(:,:,cnt)=tauy(:,:,cnt)-tauy(:,:,cnt-1); 
   end
end

cp=cpmod;
lsp=lspmod;
ssrd=ssrdmod;
strd=strdmod;
ee=eemod;
hh=hhmod;
taux=tauxmod;
tauy=tauymod;

%% soil and sea
% calculate sst
sst=zeros(ntime,1);
test=find(lsm(:,:,1)==0); % find sea in the domain
if numel(test) >0
    for i=1:ntime
       aux=esst(:,:,i);
       sst(i)=mean(aux(test));
    end    
end

% Recalculate volumetric soil moisture to the user-defined soil type
    % start with calculation psi=-psi_sat*(theta/theta_sat)^-b  (Clapp and Hornberger 1978)
    wsl1(wsl1<0) = 0;    
    wsl2(wsl2<0) = 0;
    wsl3(wsl3<0) = 0;
    wsl4(wsl4<0) = 0;
    
    psi1=-0.338*(wsl1/0.472).^(-6.04); %ECMWF parameters, see http://www.ecmwf.int/research/ifsdocs/CY37r2/IFSPart4.pdf p125
    psi2=-0.338*(wsl2/0.472).^(-6.04);
    psi3=-0.338*(wsl3/0.472).^(-6.04);
    psi4=-0.338*(wsl4/0.472).^(-6.04);
    
    wsl1=etas*(psi1/(-psis)).^(-1./bch);
    wsl2=etas*(psi2/(-psis)).^(-1./bch);
    wsl3=etas*(psi3/(-psis)).^(-1./bch);
    wsl4=etas*(psi4/(-psis)).^(-1./bch);

% calculate soil temperature and water contents
    dpth=[0.07,0.21,0.72,1.89];
    totdpth=sum(dpth);

    Ts=tsl1;
    Ts2=tsl2;
    Ts3=tsl3;
    Td=(tsl1*dpth(1)+tsl2*dpth(2)+tsl3*dpth(3)+tsl4*dpth(4)) / totdpth;

    xws=wsl1;
    xwd=(wsl1*dpth(1)+wsl2*dpth(2)+wsl3*dpth(3)+wsl4*dpth(4)) / totdpth;


    
%% extract 4x4 subdomain around clon and clat
aux=lat-clat;
aux2=find(aux>0);
lat0=lat(aux2(length(aux2))); %latitude above latc

aux=lon-clon;
aux2=find(aux<0);
lon0=lon(aux2(length(aux2))); % longitude below lonc
  
disp('lat0 and lon0')
disp(lat0)
disp(lon0)

% zoek coordinaten van dichtsbijzijne ecmwf punt
j0=find(lat==lat0);
i0=find(lon==lon0);

% coordinates of subdomain to be extracted
try
    xlon=lon(i0-1:i0+2);
    xlat=lat(j0-1:j0+2);
catch
    err = MException('InputError:CoordinatesError', ...
        'Coordinates do not match with input file.');
    throw(err)
end

disp('4x4 subdomain around city')
disp(xlon)
disp(xlat)
    
% extract variables for subdomain
xT=tt(i0-1:i0+2,j0-1:j0+2,:,:);
xq=qq(i0-1:i0+2,j0-1:j0+2,:,:);
xu=uu(i0-1:i0+2,j0-1:j0+2,:,:);
xv=vv(i0-1:i0+2,j0-1:j0+2,:,:);

Tvirt=xT.*(1.+(Rvap/Rdry-1.)*xq); % virtual temperature

xzsfc=zsfc(i0-1:i0+2,j0-1:j0+2,:);
xps=ps(i0-1:i0+2,j0-1:j0+2,:);
xgs=gs(i0-1:i0+2,j0-1:j0+2,:);
xhh=hh(i0-1:i0+2,j0-1:j0+2,:);
xee=ee(i0-1:i0+2,j0-1:j0+2,:);
xssrd=ssrd(i0-1:i0+2,j0-1:j0+2,:);
xstrd=strd(i0-1:i0+2,j0-1:j0+2,:);
xtaux=taux(i0-1:i0+2,j0-1:j0+2,:);
xtauy=tauy(i0-1:i0+2,j0-1:j0+2,:);
xlsp=lsp(i0-1:i0+2,j0-1:j0+2,:);
xcp=cp(i0-1:i0+2,j0-1:j0+2,:);
xTs=Ts(i0-1:i0+2,j0-1:j0+2,:);
xTd=Td(i0-1:i0+2,j0-1:j0+2,:);
xws=xws(i0-1:i0+2,j0-1:j0+2,:);
xwd=xwd(i0-1:i0+2,j0-1:j0+2,:);
xlsm=lsm(i0-1:i0+2,j0-1:j0+2,:);

xTs2=Ts2(i0-1:i0+2,j0-1:j0+2,:);
xTs3=Ts3(i0-1:i0+2,j0-1:j0+2,:);

%% calculate geopotential
load(modellevels_file);
levnum=modellevels(:,1);
lev_aa=modellevels(:,2);
lev_bb=modellevels(:,3);

xpp=zeros(4,4,nlevs,ntime);       % full-level pressure for subdomain
xgp=zeros(4,4,nlevs,ntime);       % full-level geopotential for subdomain

ppls=xps;   %initial pressure at NLEV+1/2 initialized with surface value
gpls=xgs;   %initial geopotential at NLEV+1/2 initialized with surface value


for kk=nlevs:-1:1
  % pressure
  aa=lev_aa(levels(kk));
  bb=lev_bb(levels(kk));
  pmin=aa+bb*xps;                 % pressure at upper interface (k-1/2)
  xpp(:,:,kk,:)=0.5*(ppls+pmin);  % full-level pressure value (level k)

  % geopotential
  alfa=1.-(pmin./(ppls-pmin)).*log(ppls./pmin);
  gmin=gpls+Rdry*squeeze(Tvirt(:,:,kk,:)).*log(ppls./pmin); % geopotential at upper interface (k-1/2)
  xgp(:,:,kk,:)=gpls+alfa.*Rdry.*squeeze(Tvirt(:,:,kk,:));   % full-level geopotential value (level k)

  % update step
  ppls=pmin; % pressure at lower interface (k+1/2)
  gpls=gmin; % geopotential at lower interface (k+1/2)    
    
end

% potential temperature
xPt=xT.*(100000./xpp).^0.287;

%% horizontal pressure gradient
xlp=log(xpp);

PGFx=zeros(4,4,nlevs,ntime);
PGFy=zeros(4,4,nlevs,ntime);

for jj=2:3
    for ii=2:3
          dgpdlon=(xgp(ii+1,jj,:,:)-xgp(ii-1,jj,:,:))/((xlon(ii+1)-xlon(ii-1))*pi/180); % zonal gradient of geopotential
          dlpdlon=(xlp(ii+1,jj,:,:)-xlp(ii-1,jj,:,:))/((xlon(ii+1)-xlon(ii-1))*pi/180); % zonal gradient of ln(pressure)
          PGFx(ii,jj,:,:)=-(dgpdlon+Rdry*Tvirt(ii,jj,:,:).*dlpdlon)/(re*cos(xlat(jj)*pi/180));

          dgpdlat=(xgp(ii,jj+1,:,:)-xgp(ii,jj-1,:,:))/((xlat(jj+1)-xlat(jj-1))*pi/180); % meridional gradient of geopotential
          dlpdlat=(xlp(ii,jj+1,:,:)-xlp(ii,jj-1,:,:))/((xlat(jj+1)-xlat(jj-1))*pi/180); % meridional gradient of ln(pressure)
          PGFy(ii,jj,:,:)=-(dgpdlat+Rdry*Tvirt(ii,jj,:,:).*dlpdlat)/re; 
          
    end
end

%% horizontal interpolation
% select the 2X2 points surrounding [clon,clat]
xlon=xlon(2:3);
xlat=xlat(2:3);
xT=xT(2:3,2:3,:,:);
xPt=xPt(2:3,2:3,:,:);
xq=xq(2:3,2:3,:,:);
xu=xu(2:3,2:3,:,:);
xv=xv(2:3,2:3,:,:);
xgp=xgp(2:3,2:3,:,:);
xgs=xgs(2:3,2:3,:);
xzsfc=xzsfc(2:3,2:3,:);
xps=xps(2:3,2:3,:);
xhh=xhh(2:3,2:3,:);
xee=xee(2:3,2:3,:);
xssrd=xssrd(2:3,2:3,:);
xstrd=xstrd(2:3,2:3,:);
xtaux=xtaux(2:3,2:3,:);
xtauy=xtauy(2:3,2:3,:);
xlsp=xlsp(2:3,2:3,:);
xcp=xcp(2:3,2:3,:);
xTs=xTs(2:3,2:3,:);
xTd=xTd(2:3,2:3,:);
xws=xws(2:3,2:3,:);
xwd=xwd(2:3,2:3,:);
xlsm=xlsm(2:3,2:3,:);
xTs2=xTs2(2:3,2:3,:);
xTs3=xTs3(2:3,2:3,:);
PGFx=PGFx(2:3,2:3,:,:);
PGFy=PGFy(2:3,2:3,:,:);

% bilinear interpolation
ic=(clon-xlon(1))/(xlon(2)-xlon(1));
jc=(clat-xlat(1))/(xlat(2)-xlat(1));

cT=xT(1,1,:,:)*(1.-ic)*(1.-jc)+xT(2,1,:,:)*ic*(1.-jc)+xT(1,2,:,:)*(1.-ic)*jc+xT(2,2,:,:)*ic*jc;
cPt=xPt(1,1,:,:)*(1.-ic)*(1.-jc)+xPt(2,1,:,:)*ic*(1.-jc)+xPt(1,2,:,:)*(1.-ic)*jc+xPt(2,2,:,:)*ic*jc;
cq=xq(1,1,:,:)*(1.-ic)*(1.-jc)+xq(2,1,:,:)*ic*(1.-jc)+xq(1,2,:,:)*(1.-ic)*jc+xq(2,2,:,:)*ic*jc;
cu=xu(1,1,:,:)*(1.-ic)*(1.-jc)+xu(2,1,:,:)*ic*(1.-jc)+xu(1,2,:,:)*(1.-ic)*jc+xu(2,2,:,:)*ic*jc;
cv=xv(1,1,:,:)*(1.-ic)*(1.-jc)+xv(2,1,:,:)*ic*(1.-jc)+xv(1,2,:,:)*(1.-ic)*jc+xv(2,2,:,:)*ic*jc;
cgp=xgp(1,1,:,:)*(1.-ic)*(1.-jc)+xgp(2,1,:,:)*ic*(1.-jc)+xgp(1,2,:,:)*(1.-ic)*jc+xgp(2,2,:,:)*ic*jc;
cgs=xgs(1,1,:)*(1.-ic)*(1.-jc)+xgs(2,1,:)*ic*(1.-jc)+xgs(1,2,:)*(1.-ic)*jc+xgs(2,2,:)*ic*jc;
czsfc=xzsfc(1,1,:)*(1.-ic)*(1.-jc)+xzsfc(2,1,:)*ic*(1.-jc)+xzsfc(1,2,:)*(1.-ic)*jc+xzsfc(2,2,:)*ic*jc;
cps=xps(1,1,:)*(1.-ic)*(1.-jc)+xps(2,1,:)*ic*(1.-jc)+xps(1,2,:)*(1.-ic)*jc+xps(2,2,:)*ic*jc;
chh=xhh(1,1,:)*(1.-ic)*(1.-jc)+xhh(2,1,:)*ic*(1.-jc)+xhh(1,2,:)*(1.-ic)*jc+xhh(2,2,:)*ic*jc;
cee=xee(1,1,:)*(1.-ic)*(1.-jc)+xee(2,1,:)*ic*(1.-jc)+xee(1,2,:)*(1.-ic)*jc+xee(2,2,:)*ic*jc;
cssrd=xssrd(1,1,:)*(1.-ic)*(1.-jc)+xssrd(2,1,:)*ic*(1.-jc)+xssrd(1,2,:)*(1.-ic)*jc+xssrd(2,2,:)*ic*jc;
cstrd=xstrd(1,1,:)*(1.-ic)*(1.-jc)+xstrd(2,1,:)*ic*(1.-jc)+xstrd(1,2,:)*(1.-ic)*jc+xstrd(2,2,:)*ic*jc;
ctaux=xtaux(1,1,:)*(1.-ic)*(1.-jc)+xtaux(2,1,:)*ic*(1.-jc)+xtaux(1,2,:)*(1.-ic)*jc+xtaux(2,2,:)*ic*jc;
ctauy=xtauy(1,1,:)*(1.-ic)*(1.-jc)+xtauy(2,1,:)*ic*(1.-jc)+xtauy(1,2,:)*(1.-ic)*jc+xtauy(2,2,:)*ic*jc;
clsp=xlsp(1,1,:)*(1.-ic)*(1.-jc)+xlsp(2,1,:)*ic*(1.-jc)+xlsp(1,2,:)*(1.-ic)*jc+xlsp(2,2,:)*ic*jc;
ccp=xcp(1,1,:)*(1.-ic)*(1.-jc)+xcp(2,1,:)*ic*(1.-jc)+xcp(1,2,:)*(1.-ic)*jc+xcp(2,2,:)*ic*jc;
cTs=xTs(1,1,:)*(1.-ic)*(1.-jc)+xTs(2,1,:)*ic*(1.-jc)+xTs(1,2,:)*(1.-ic)*jc+xTs(2,2,:)*ic*jc;
cTd=xTd(1,1,:)*(1.-ic)*(1.-jc)+xTd(2,1,:)*ic*(1.-jc)+xTd(1,2,:)*(1.-ic)*jc+xTd(2,2,:)*ic*jc;
clsm=xlsm(1,1,:)*(1.-ic)*(1.-jc)+xlsm(2,1,:)*ic*(1.-jc)+xlsm(1,2,:)*(1.-ic)*jc+xlsm(2,2,:)*ic*jc;
cTs2=xTs2(1,1,:)*(1.-ic)*(1.-jc)+xTs2(2,1,:)*ic*(1.-jc)+xTs2(1,2,:)*(1.-ic)*jc+xTs2(2,2,:)*ic*jc;
cTs3=xTs3(1,1,:)*(1.-ic)*(1.-jc)+xTs3(2,1,:)*ic*(1.-jc)+xTs3(1,2,:)*(1.-ic)*jc+xTs3(2,2,:)*ic*jc;
cPGFx=PGFx(1,1,:,:)*(1.-ic)*(1.-jc)+PGFx(2,1,:,:)*ic*(1.-jc)+PGFx(1,2,:,:)*(1.-ic)*jc+PGFx(2,2,:,:)*ic*jc;
cPGFy=PGFy(1,1,:,:)*(1.-ic)*(1.-jc)+PGFy(2,1,:,:)*ic*(1.-jc)+PGFy(1,2,:,:)*(1.-ic)*jc+PGFy(2,2,:,:)*ic*jc;

test=find(xws(:,:,1)>0.03);
if numel(test)>0
   for i=1:ntime
       aux=xws(:,:,i);
       cws(i)=mean(aux(test));
       aux2=xwd(:,:,i);
       cwd(i)=mean(aux2(test));
   end 
    
end


cT=squeeze(cT);
cPt=squeeze(cPt);
cq=squeeze(cq);
cu=squeeze(cu);
cv=squeeze(cv);
cgp=squeeze(cgp);
cgs=squeeze(cgs);
czsfc=squeeze(czsfc);
cps=squeeze(cps);
chh=squeeze(chh);
cee=squeeze(cee);
cssrd=squeeze(cssrd);
cstrd=squeeze(cstrd);
ctaux=squeeze(ctaux);
ctauy=squeeze(ctauy);
clsp=squeeze(clsp);
ccp=squeeze(ccp);
cTs=squeeze(cTs);
cTd=squeeze(cTd);
clsm=squeeze(clsm);
cws=squeeze(cws);
cwd=squeeze(cwd);
cTs2=squeeze(cTs2);
cTs3=squeeze(cTs3);
cPGFx=squeeze(cPGFx);
cPGFy=squeeze(cPGFy);



%% vertical interpolation
nout=length(zout);
Tout1=zeros(nout,ntime);
Ptout1=zeros(nout,ntime);
qout1=zeros(nout,ntime);
uout1=zeros(nout,ntime);
vout1=zeros(nout,ntime);

PGFxout1=zeros(nout,ntime);
PGFyout1=zeros(nout,ntime);


for ii=1:ntime

  cz=(cgp(:,ii)-cgs(ii))/grav;  
   if min(cz) > zout(1)
       verschil = cz(length(cz)) -zout(1);
       cz(length(cz))=zout(1);
       if verschil > 1
           cz
           zout
           disp(verschil)
           err = MException('InputError:ECMWFlevels', ...
            'levels in input data file do not extend low enough');
            throw(err)            
       end
       
   end

  if max(cz) < max(zout) 
    disp 'levels in input data file do not extend high enough'
    disp(ii)
    disp(max(cz))
    disp(max(zout))
    disp 'either modify input data, or select lower UCM model top'
    err = MException('InputError:InputECMWF', ...
        'levels in input data file do not extend high enough');
    throw(err)    
  end
  
  Tout1(:,ii)=interp1(cz,cT(:,ii),zout);
  Ptout1(:,ii)=interp1(cz,cPt(:,ii),zout);
  qout1(:,ii)=interp1(cz,cq(:,ii),zout);
  uout1(:,ii)=interp1(cz,cu(:,ii),zout);
  vout1(:,ii)=interp1(cz,cv(:,ii),zout);

  PGFxout1(:,ii)=interp1(cz,cPGFx(:,ii),zout);
  PGFyout1(:,ii)=interp1(cz,cPGFy(:,ii),zout);

end

%% temporal interpolation
% resample from 3-hourly to new resolution

tmet_day=tmet/(3600*24);
time_output=matlabtime(1):tmet_day:matlabtime(length(matlabtime));
matlabtime=unique(matlabtime);
xx=length(time_output);

Tout=zeros(nout,xx);
Ptout=zeros(nout,xx);
qout=zeros(nout,xx);
uout=zeros(nout,xx);
vout=zeros(nout,xx);
PGFxout=zeros(nout,xx);
PGFyout=zeros(nout,xx);

for ii=1:nout
  Tout(ii,:)=interp1(matlabtime,Tout1(ii,:),time_output);
  Ptout(ii,:)=interp1(matlabtime,Ptout1(ii,:),time_output);
  qout(ii,:)=interp1(matlabtime,qout1(ii,:),time_output);
  uout(ii,:)=interp1(matlabtime,uout1(ii,:),time_output);
  vout(ii,:)=interp1(matlabtime,vout1(ii,:),time_output);

  PGFxout(ii,:)=interp1(matlabtime,PGFxout1(ii,:),time_output);
  PGFyout(ii,:)=interp1(matlabtime,PGFyout1(ii,:),time_output);

end

ps = interp1(matlabtime,cps,time_output);
hh = interp1(matlabtime,chh,time_output);
ee =interp1(matlabtime,cee,time_output);
ssrd = interp1(matlabtime,cssrd,time_output);
strd = interp1(matlabtime,cstrd,time_output);
taux = interp1(matlabtime,ctaux,time_output);
tauy = interp1(matlabtime,ctauy,time_output);
lsp = interp1(matlabtime,clsp,time_output);
cp = interp1(matlabtime,ccp,time_output);
Ts = interp1(matlabtime,cTs,time_output);
Td = interp1(matlabtime,cTd,time_output);
Ws = interp1(matlabtime,cws,time_output);
Wd = interp1(matlabtime,cwd,time_output);
lsm = interp1(matlabtime,clsm,time_output);
sst = interp1(matlabtime,sst,time_output);

Ts2 = interp1(matlabtime,cTs2,time_output);
Ts3 = interp1(matlabtime,cTs3,time_output);

%% finishing
zsfc_cen=czsfc(1);
disp 'surface height is'
disp(zsfc_cen)

rain=(lsp+cp)/(10800./tmet);   %due to interpolation from 3 hourly to user-defined resolution
stress=sqrt(taux.^2.+tauy.^2.);
tvirt=Tout(1,:).*(1.+0.61.*qout(1,:));    % virtual temperature
rho=ps./(287.*tvirt);       % air density
ustr=sqrt(stress./rho);    % friction velocity

year=str2num(datestr(time_output,'yyyy'));
month=str2num(datestr(time_output,'mm'));
day=str2num(datestr(time_output,'dd'));
hour=str2num(datestr(time_output,'HH'));
minute=str2num(datestr(time_output,'MM'));


%% visualisation
% enkele tijdsreeksen
aantal_reeksen=3;

f=figure('visible','off');
hold on
set(f,'Position',[0,0,1000,1000]);
subplot(aantal_reeksen,1,1);
plot(time_output,hh);
datetick('x','mm-YYYY','keeplimits')
ylabel('H')
subplot(aantal_reeksen,1,2);
plot(time_output,ee);
datetick('x','mm-YYYY','keeplimits')
ylabel('LE')
subplot(aantal_reeksen,1,3);
plot(time_output,ssrd);
datetick('x','mm-YYYY','keeplimits')
ylabel('Rs')
print(f,'-dpng',strcat(logfolder,pretext,'ECMWF_timeseries_matlab1.png'));
close(f)

f=figure('visible','off');
hold on
set(f,'Position',[0,0,1000,1000]);
subplot(aantal_reeksen,1,1);
plot(time_output,strd);
datetick('x','mm-YYYY','keeplimits')
ylabel('Rl')
subplot(aantal_reeksen,1,2);
plot(time_output,ustr);
datetick('x','mm-YYYY','keeplimits')
ylabel('u*')
subplot(aantal_reeksen,1,3);
plot(time_output,ps/10);
datetick('x','mm-YYYY','keeplimits')
ylabel('ps')
print(f,'-dpng',strcat(logfolder,pretext,'ECMWF_timeseries_matlab2.png'));
close(f)


f=figure('visible','off');
hold on
set(f,'Position',[0,0,1000,1000]);
subplot(aantal_reeksen,1,1);
hold on
plot(time_output,Ts);
plot(time_output,Td,'k');
plot(time_output,Ts3,'r');
datetick('x','mm-YYYY','keeplimits')
ylabel('Ts')
subplot(aantal_reeksen,1,2);
hold on
plot(time_output,Ws);
plot(time_output,Wd,'r');
datetick('x','mm-YYYY','keeplimits')
ylabel('Ws')
subplot(aantal_reeksen,1,3);
plot(time_output,rain);
datetick('x','mm-YYYY','keeplimits')
ylabel('rain')
print(f,'-dpng',strcat(logfolder,pretext,'ECMWF_timeseries_matlab3.png'));
close(f)

f=figure('visible','off');
hold on
set(f,'Position',[0,0,1000,1000]);
subplot(3,1,1);
plot(time_output,sst);
ylabel('SST')
datetick('x','mm-YYYY','keeplimits')
subplot(3,1,2);
plot(time_output,Tout(1,:));
ylabel('Tsurface')
datetick('x','mm-YYYY','keeplimits')
subplot(3,1,3);
surfacewind=sqrt(uout(1,:).^2+vout(1,:).^2);
plot(time_output,surfacewind(1,:));
ylabel('Wind surface')
datetick('x','mm-YYYY','keeplimits')
print(f,'-dpng',strcat(logfolder,pretext,'ECMWF_timeseries_matlab4.png'));
close(f)


% enkele verticale temperatuursprofielen
f=figure('visible','off');
hold on
set(f,'Position',[0,0,1000,1000]);


for i=1:8
    subplot(2,4,i)
    moment=28+(i-1)*36;
    plot(Tout(:,moment),zout);
    axis([min(Tout(:,moment)) max(Tout(:,moment)) 0 2500]);
    title([num2str(day(moment)),'/',num2str(month(moment)),'/',num2str(year(moment)),' T',num2str(hour(moment)),'h']);  
end
print(f,'-dpng','-r300',strcat(logfolder,pretext,'ECMWF_vertical_profiles_matlab.png'));
close(f)



%% write data
    id = netcdf.create(outputfile,'CLOBBER');

    idtimedim = netcdf.defDim(id,'time',length(time_output));
    idleveldim = netcdf.defDim(id,'nz',nout);
       
    idtimevar = netcdf.defVar(id,'matlabtime','double',idtimedim);
    idlevelvar = netcdf.defVar(id,'level','long',idleveldim);    
    
    vid0 = netcdf.defVar(id, 'YEAR','long',idtimedim);
    vid1 = netcdf.defVar(id, 'MONTH', 'long',idtimedim);
    vid2 = netcdf.defVar(id, 'DAY', 'long',idtimedim);
    vid3 = netcdf.defVar(id, 'HOUR', 'long',idtimedim);
    vid4 = netcdf.defVar(id, 'MINUTE', 'long',idtimedim);
    vid5 = netcdf.defVar(id, 'LATC','float',[]);
    vid6 = netcdf.defVar(id, 'LONC', 'float',[]);
    vid7 = netcdf.defVar(id, 'Altitude', 'float',[]);
    vid8 = netcdf.defVar(id, 'Z', 'float',idleveldim);
    vid9 = netcdf.defVar(id, 'H',  'float',idtimedim);
    vid10 = netcdf.defVar(id, 'LE',  'float',idtimedim);
    vid11 = netcdf.defVar(id, 'PS', 'float',idtimedim);
    vid12 = netcdf.defVar(id, 'Ts', 'float',idtimedim);
    vid13 = netcdf.defVar(id, 'Td', 'float',idtimedim);
    vid14 = netcdf.defVar(id, 'Tw', 'float',idtimedim);
    vid15 = netcdf.defVar(id, 'SST', 'float',idtimedim);
    vid16 = netcdf.defVar(id, 'Ws', 'float',idtimedim);
    vid17 = netcdf.defVar(id, 'Wd', 'float',idtimedim);
    vid18 = netcdf.defVar(id, 'Rain', 'float',idtimedim);
    vid19 = netcdf.defVar(id, 'Rs', 'float',idtimedim);
    vid20 = netcdf.defVar(id, 'Rl','float',idtimedim);
    vid21 = netcdf.defVar(id, 'U', 'float',[idleveldim,idtimedim]);
    vid22 = netcdf.defVar(id, 'V', 'float',[idleveldim,idtimedim]);
    vid23 = netcdf.defVar(id, 'T', 'float',[idleveldim,idtimedim]);
    vid24 = netcdf.defVar(id, 'Pt', 'float',[idleveldim,idtimedim]);
    vid25 = netcdf.defVar(id, 'Q', 'float',[idleveldim,idtimedim]);
    vid26 = netcdf.defVar(id, 'PGFx', 'float',[idleveldim,idtimedim]);
    vid27 = netcdf.defVar(id, 'PGFy', 'float',[idleveldim,idtimedim]);
    
    netcdf.putAtt(id,vid0,'long_name', 'Year');
    netcdf.putAtt(id,vid1,'long_name', 'Month');
    netcdf.putAtt(id,vid2,'long_name', 'Day');
    netcdf.putAtt(id,vid3,'long_name', 'Hour');
    netcdf.putAtt(id,vid4,'long_name', 'Minute');
    netcdf.putAtt(id,vid5,'long_name', 'Latitude');
    netcdf.putAtt(id,vid5,'units','degrees');
    netcdf.putAtt(id,vid6,'long_name', 'Longitude');
    netcdf.putAtt(id,vid6,'units','degrees');
    netcdf.putAtt(id,vid7,'long_name', 'Altitude');
    netcdf.putAtt(id,vid7,'units','meters');
    netcdf.putAtt(id,vid8,'long_name', 'Height of the model levels');
    netcdf.putAtt(id,vid8,'units','meters');
    netcdf.putAtt(id,vid9,'long_name', 'Sensible heat flux');
    netcdf.putAtt(id,vid9,'units','W/m2');
    netcdf.putAtt(id,vid10,'long_name', 'Latent heat flux');
    netcdf.putAtt(id,vid10,'units','W/m2');
    netcdf.putAtt(id,vid11,'long_name', 'Surface pressure');
    netcdf.putAtt(id,vid11,'units','hPa');
    netcdf.putAtt(id,vid12,'long_name', 'Skin temperature');
    netcdf.putAtt(id,vid12,'units','K');
    netcdf.putAtt(id,vid13,'long_name', 'Deep soil temperature');
    netcdf.putAtt(id,vid13,'units','K');
    netcdf.putAtt(id,vid14,'long_name', 'Temperature for water bodies');
    netcdf.putAtt(id,vid14,'units','K');
    netcdf.putAtt(id,vid15,'long_name', 'Sea surface temperature');
    netcdf.putAtt(id,vid15,'units','K');
    netcdf.putAtt(id,vid16,'long_name', 'Surface wetness');
    netcdf.putAtt(id,vid16,'units','percent');
    netcdf.putAtt(id,vid17,'long_name', 'Deep soil wetness');
    netcdf.putAtt(id,vid17,'units','percent');
    netcdf.putAtt(id,vid18,'long_name', 'Rainfall');
    netcdf.putAtt(id,vid18,'units','mm');
    netcdf.putAtt(id,vid19,'long_name', 'Incoming solar radiation');
    netcdf.putAtt(id,vid19,'units','W/m2');
    netcdf.putAtt(id,vid20,'long_name', 'Incoming longwave radiation');
    netcdf.putAtt(id,vid20,'units','W/m2');
    netcdf.putAtt(id,vid21,'long_name', 'Horizontal wind speed in x');
    netcdf.putAtt(id,vid21,'units','m/s');
    netcdf.putAtt(id,vid22,'long_name', 'Horizontal wind speed in y');
    netcdf.putAtt(id,vid22,'units','m/s');
    netcdf.putAtt(id,vid23,'long_name', 'Temperature');
    netcdf.putAtt(id,vid23,'units','K');
    netcdf.putAtt(id,vid24,'long_name', 'Potential temperature');
    netcdf.putAtt(id,vid24,'units','K');
    netcdf.putAtt(id,vid25,'long_name', 'Specific humidity');
    netcdf.putAtt(id,vid25,'units','g/kg');
    netcdf.putAtt(id,vid26,'long_name', 'Pressure gradient force in x');
    netcdf.putAtt(id,vid26,'units','m/s2');
    netcdf.putAtt(id,vid27,'long_name', 'Pressure gradient force in y');
    netcdf.putAtt(id,vid27,'units','m/s2');

    netcdf.endDef(id);
    
    netcdf.putVar(id,idtimevar,time_output);
    netcdf.putVar(id,idlevelvar,1:length(zout));
    
    netcdf.putVar(id,vid0,year);
    netcdf.putVar(id,vid1,month);
    netcdf.putVar(id,vid2,day);
    netcdf.putVar(id,vid3,hour);
    netcdf.putVar(id,vid4,minute);
    netcdf.putVar(id,vid5,clat);
    netcdf.putVar(id,vid6,clon);
    netcdf.putVar(id,vid7,zsfc_cen);
    netcdf.putVar(id,vid8,zout);
    netcdf.putVar(id,vid9,hh);
    netcdf.putVar(id,vid10,ee);
    netcdf.putVar(id,vid11,ps/100.); %in hPa
    netcdf.putVar(id,vid12,Ts);
    netcdf.putVar(id,vid13,Td);
    netcdf.putVar(id,vid14,Ts3);
    netcdf.putVar(id,vid15,sst);
    netcdf.putVar(id,vid16,Ws);
    netcdf.putVar(id,vid17,Wd);
    netcdf.putVar(id,vid18,rain);
    netcdf.putVar(id,vid19,ssrd);
    netcdf.putVar(id,vid20,strd);
    netcdf.putVar(id,vid21,uout);
    netcdf.putVar(id,vid22,vout);
    netcdf.putVar(id,vid23,Tout);
    netcdf.putVar(id,vid24,Ptout);
    netcdf.putVar(id,vid25,qout*1000.); %in g/kg
    netcdf.putVar(id,vid26,PGFxout);
    netcdf.putVar(id,vid27,PGFyout);
    
    netcdf.close(id);
        
disp('Done')

end

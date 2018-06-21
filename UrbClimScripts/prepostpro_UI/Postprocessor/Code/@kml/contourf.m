function target = contourf(this,long,lat,alt,levels,colordata,varargin)
%KML.CONTOURF(long,lat,alt) Create a filled contour of alt in a grid defined by long and lat. 
%   Similar to built-in contour function
%
%   17.05.2013 - Added suggestion by Keith Epstein to allow changing line width and color
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $


    target = struct('type','','id','');
    
    p = inputParser;
    
    nlat = numel(lat);

    p.addRequired('lat',  @(a)isnumeric(a) && ~isempty(a));
    p.addRequired('long', @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('alt',  @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('levels',  @(a)isnumeric(a) && ~isempty(a));
    
    
    
    p.addParamValue('name','kml_contourf',@ischar);
    p.addParamValue('id',kml.getTempID('kml_contourf'),@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('transparency',1,@isnumeric);
    p.addParamValue('lineWidth',1,@isnumeric);
    p.addParamValue('lineColor','',@(x)ischar(x)||isempty(x));
    
    p.addParamValue('noFolder',false,@islogical)
    
    p.addParamValue('showText',false,@islogical)
    p.addParamValue('levelStep',1,@isnumeric)
    p.addParamValue('labelSpacing',inf,@isnumeric)

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,alt,levels,varargin{:});
    
    arg = p.Results;
    
    if arg.noFolder
        f = this;
    else
        f = this.createFolder(arg.name);
    end
            
    
    [m,n] = size(alt);
    
    alt = [NaN(1, n+2); NaN(m,1) alt NaN(m,1); NaN(1, n+2)];
    lat = [2.*lat(1,[1 1:end end]) - lat(2,[1 1:end end]); ...
           (2.*lat(:,1)-lat(:,2)), lat, (2.*lat(:,end)-lat(:,end-1)); ...
           2.*lat(end,[1 1:end end]) - lat(end-1,[1 1:end end])];
       
    long = long.';
    long = [2.*long(1,[1 1:end end]) - long(2,[1 1:end end]); ...
           (2.*long(:,1)-long(:,2)), long, (2.*long(:,end)-long(:,end-1)); ...
           2.*long(end,[1 1:end end]) - long(end-1,[1 1:end end])];
       
    long = long.';
    

    
    alt(isnan(alt)) = min(levels) - 1e4.*(max(levels) - min(levels));    
    
    alt(m+2,m+2) = max(levels) + 1e4.*(max(levels) - min(levels));  
    
    c = contours(long,lat,alt,levels);
    
    
    i = 1; k=1;
    while true
        l = c(1,i);
        n = c(2,i);
        C(k).long = c(1,(i+1):(i+n));
        C(k).lat = c(2,(i+1):(i+n));
        C(k).l = l;
        if n<2
            C(k).area = 0;
        else
            C(k).area = sum(diff(C(k).lat) .* (C(k).long(1 : n - 1) + C(k).long(2 : n)) / 2);
        end

        
        i = i + n + 1;
        k = k + 1;
        if i > size(c,2)
            break;
        end
    end
    

    levels = unique([C(:).l]);
    minLevel = min(levels);
    shownLevels = levels(1:arg.levelStep:numel(levels));
    areas = -unique(-abs([C(:).area]));
    [~,ixMA] = max(abs([C(:).area]));
    maxArea = C(ixMA).area;
    
    ncolors = 100;
    cmap = colordata;
    ncolors=length(cmap);
    levelsCMAP = linspace(min(levels),max(levels),ncolors);
    middleLevel = levels(ceil(numel(levels)/2));
    
    
    % Don't fill contours below the minimum level
    showMinLevel=any(levels <= min(levels));
    
    % Don't fill contours above the maximum level
    showMaxLevel=any(levels >= max(levels));
    
    for a = areas;
        for i = 1:numel(C)
            if abs(C(i).area) == a
                clev = C(i).l;
                %iC = find(levels==clev);
%                 if clev <= middleLevel
%                     clev = clev - 1;
%                 end

                if (sign(C(i).area) ~=sign(maxArea)),
                    kk=find(levels==clev);
                    kk0 = 1 + sum(levels<=min(levels)) * (~showMinLevel);
                    if (kk > kk0)
                        clev=levels(kk-1);    % in valley, use color for lower level
                    elseif (kk == kk0)
                        break;
                    else
                        break;
                    end
                end

                if numel(levels) > 1
                    iC = floor(interp1(levelsCMAP,linspace(0,ncolors-1,numel(levelsCMAP)),clev,'linear',0));
                else
                    iC = 1;
                end

                color = cmap(iC+1 ,:);
                
                colorHex = kml.color2kmlHex([color arg.transparency]);
                
                % debug: to view while plotting
                % fig(12345)
                % hold on
                % patch(C(i).long,C(i).lat,color)
                
                
                if arg.lineWidth > 0
                    lineColor = '00000000';
                else
                    if isempty(arg.lineColor)
                        lineColor = ['FF' colorHex(3:end)];
                    else
                        lineColor = 'FF000000';
                    end
                end
                
                target(end+1) = f.poly3(C(i).long,C(i).lat,zeros(size(C(i).lat)), 'polyColor', colorHex, ...
                                           'lineColor',lineColor,...
                                           'lineWidth',arg.lineWidth, ...
                                           'altitudeMode','clampToGround', ...
                                           'visibility',arg.visibility, ...
                                           'name',sprintf('Level %g',C(i).l), ...
                                           'timeStamp', arg.timeStamp , ...
                                           'timeSpanBegin', arg.timeSpanBegin , ...
                                           'timeSpanEnd', arg.timeSpanEnd, ...      
                                           'id',[arg.id '_poly_' num2str(i)] ...  
                                           );
                if arg.showText && ismember(C(i).l,shownLevels)
                    N = numel(C(i).lat);
                    if ~isfinite(arg.labelSpacing)
                        R = randi(N,1);
                        latTxt  = interp1(1:N,C(i).lat,R);
                        longTxt = interp1(1:N,C(i).long,R);
                        target(end+1) = f.text(longTxt,latTxt,0,sprintf('%g',C(i).l), 'id',[arg.id '_text_' num2str(i)]);
                    else
                       dist = kml.ll2dist(C(i).long(1:end-1),C(i).lat(1:end-1), ...
                                          C(i).long(2:end),C(i).lat(2:end));
                       dist = [0 cumsum(dist)];
                       for d = 1:arg.labelSpacing:dist(end);
                           latTxt  = interp1(dist,C(i).lat,d);
                           longTxt = interp1(dist,C(i).long,d);
                           target(end+1) = f.text(longTxt,latTxt,0,sprintf('%g',C(i).l), 'id',[arg.id '_text_' num2str(i) '_' num2str(d)]);
                       end
                    end
                end
            end
        end
    end    
    target(1) = []; %remove the empty initial field
end

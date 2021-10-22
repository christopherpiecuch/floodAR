%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to download tidal datums from NOAA NOS
% Original code written by Dr. Steven J. Lentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [datum]=noaaDatums(station)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [datum]=noaaDatums(station)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downloads datums from NOAA at tide gauge 
% given by "station" ID (e.g., WHOI is 8447930) 
% all values are in m relative to 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first get metadata for station
clear URL filename status X
URL=['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=20130808%2015:00&end_date=20130808%2015:06&station=',num2str(station),'&product=water_temperature&units=english&time_zone=gmt&application=ports_screen&format=xml'];
filename='data.txt';
[~,status]=urlwrite(URL, filename);
if (status==1)
    % load metadata
    X=importdata(filename);
    X=char(X(3));
    
    % now find where the metadata are
    % they're all enclosed within quotes
    for n=1:numel(X)
        if strcmp(X(n),'"')
            N(n)=1;
        else
            N(n)=0;
        end
    end
    N=find(N==1);
    N(1:2:end)=N((1:2:end))+1;
    N(2:2:end)=N((2:2:end))-1;
    
    % grab id, name, lat, lon
    datum.id=str2num(X(N(1):N(2)));
    datum.name=(X(N(3):N(4)));
    datum.lat=str2num(X(N(5):N(6)));
    datum.lon=str2num(X(N(7):N(8)));
else
    datum.id=[];
    datum.name=[];
    datum.lat=[];
    datum.lon=[];
end
 
URL=['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?station=',num2str(station),'&product=datums&units=metric&pplication=ports_screen&format=csv'];
filename='data.txt';
[~,status]=urlwrite(URL, filename);
if (status==1)
    X=importdata(filename);
    for nn=1:numel(X.data)
        eval(['datum.',char(X.textdata(nn+1,1)),'=X.data(',num2str(nn),');'])
    end
    datum.units='m STND';
else
    datum.units=[];
end


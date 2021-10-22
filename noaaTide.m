%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to download hourly tidal predictions from NOAA NOS
% Original code written by Dr. Steven J. Lentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sl,md]=noaaTide(station,start_yr,end_yr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [sl,md]=noaaSealevel(station,start_yr,end_yr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downloads hourly tidal prediction from NOAA at tide gauge 
% given by "station" ID (e.g., WHOI is 8447930) between start_yr:end_yr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize date numbers
md1=[]; md2=[]; dt=[];
dt=1/24;
md1=datenum(start_yr,1,1,0,0,0);
md2=datenum(end_yr,12,31,23,0,0);
md=md1:dt:md2;
clear md1 md2

% time parameters
iyr=start_yr:1:end_yr;
nyr=length(iyr);

% initialize tidal prediction loop variables
mdi=[];sli=[];

% loop over years
for k=1:nyr;
    yr=num2str(iyr(k));
    URL=['https://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=',...
      yr,'0101&end_date=',yr,'1231&datum=MHHW&station=',int2str(station),'&time_zone=GMT&units=metric&format=csv&interval=60'];
    filename='data.txt';
    [~,status]=urlwrite(URL, filename);
    if (status==1)
        X=importdata(filename);
        if numel(X)==1
            md1=datenum(X.textdata(2:end,1));
            sl1=X.data(:,1);
            mdi=[mdi;md1];
            sli=[sli;sl1];
        else
            X=char(X);
            X(1,:)=[];
            md1=datenum(X(:,1:16));
            sl1=str2num(X(:,18:23));
            mdi=[mdi;md1];
            sli=[sli;sl1];
        end
    end
end
delete(filename);
%put on a constant time base

if ~isempty(sli)
    %dt=1/24;
    %md=(mdi(1):dt:mdi(end)).';
    sl=nan*md;
    k=round((mdi-md(1))./dt)+1;
    sl(k)=sli;
    readme='md matlab datenumber GMT; sl sea level meters';
else
    sl=nan;
    md=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity code; used to produce different lines in Figure 3a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

% loop through different criteria for determining whether an AR is
% nearby a tide gauge (IVT 250 or 500 kg/m/s; nearest 1 or 2
% grid cells)
for nnn=1:4

% set thresholds
timeThreshold=24; % 24 hours
distThreshold=1^(mod(nnn,2)==1)+2^(mod(nnn,2)==0)-1;  % 1 or 2 nearest grid cell(s)
trspThreshold=500^(nnn<=2)+250^(nnn>2)-1; % 250 or 500 kg/m/s

% load sio r1 catalog
load('SIO_R1_1948-2017_Comprehensive.mat')

% delete storms with peak ivt < trspThreshold
clear nn ii
for nn=1:numel(unique(Number))
 ii=[]; ii=find(Number==nn);
 MM(nn)=max(IVT(ii));
end
clear nn ii
ii=[]; ii=find(MM<trspThreshold); % delete
jj=[]; jj=find(ismember(Number,ii));
Day(jj)=[];
Hour(jj)=[];
IVT(jj)=[];
IWT(jj)=[];
Latitude(jj)=[];
Longitude(jj)=[];
Month(jj)=[];
Number(jj)=[];
Uwind(jj)=[];
Vwind(jj)=[];
Year(jj)=[];
clear ii jj

jj=[]; jj=find(str2num(datestr(datenum(Year,Month,Day,Hour,0,0),10))<1980|str2num(datestr(datenum(Year,Month,Day,Hour,0,0),10))>2016);
Day(jj)=[];
Hour(jj)=[];
IVT(jj)=[];
IWT(jj)=[];
Latitude(jj)=[];
Longitude(jj)=[];
Month(jj)=[];
Number(jj)=[];
Uwind(jj)=[];
Vwind(jj)=[];
Year(jj)=[];
clear ii jj

Longitude=Longitude-360; % degrees west

% at this point you have unique storms that meet the transport/ivt threshold
% that happened after 1980

% redefine event numbers
uNum=unique(Number);
for nn=1:numel(uNum)
 Number(find(Number==uNum(nn)))=nn;
end

clearvars -except Day Hour IVT IWT Latitude Longitude Month Number Uwind Vwind Year *Threshold

% define tide gauge identifiers
load fileID.mat

uCoords=unique(Latitude+i*Longitude);
uLat=real(uCoords);
uLon=imag(uCoords);
catalogueTime=datenum(Year,Month,Day,Hour,0,0);
hourPerDay=24;
avDa=10*2; % 20-day filtering period

for nn=1:numel(ID), disp(['Tide Gauge ',num2str(nn)])
 clearvars -except Day Hour IVT IWT Latitude Longitude Month Number Uwind Vwind Year *Threshold uCoords uLat uLon catalogueTime hourPerDay nn ID avDa

 % load data and compute tide plus mean sea level component
 load(['noaa_tidegauge_',num2str(ID(nn)),'.mat'])
 mt=td+movmedian(sl-td,hourPerDay*avDa,'omitnan','Endpoints','shrink');
 mt(find(isnan(sl)))=nan;

 disp([num2str(nn),': ',datum.name])
 sl(find(str2num(datestr(dn,10))<1980|str2num(datestr(dn,10))>2016))=[];
 td(find(str2num(datestr(dn,10))<1980|str2num(datestr(dn,10))>2016))=[];
 mt(find(str2num(datestr(dn,10))<1980|str2num(datestr(dn,10))>2016))=[];
 dn(find(str2num(datestr(dn,10))<1980|str2num(datestr(dn,10))>2016))=[];

 % identify closest catalog grid cell
 dd=[]; dd=distance(Latitude,Longitude,datum.lat,datum.lon,6371);
 sd=sort(unique(dd));
 jj=[]; jj=find(dd<=sd(distThreshold)); % take from 2 nearest grid cells

 % identify minor high-tide floods and whether they were expected from the tides+msl
 % derive flood threshold after Sweet et al. (2018)
 derFloThr=0.04*datum.GT+0.50;
 floodHour=[]; floodHour=dn(find(sl>=derFloThr))';
 tidalHour=[]; tidalHour=dn(find(mt>=derFloThr))';
 floodTide=zeros(size(floodHour));
 for ll=1:numel(floodHour);
  floodTide(ll)=floodTide(ll)+double(ismember(floodHour(ll),tidalHour));
 end 

 floodHourNearestStorm=nan*floodHour;
 % identify how close in time the nearest atmospheric river was
 for mm=1:numel(floodHour)
  floodHourNearestStorm(mm)=hourPerDay*min(abs(floodHour(mm)-catalogueTime(jj)));
 end

 % identify all atmospheric rivers
 stormHour=[]; stormHour=catalogueTime(jj);
 % identify how close in time was the nearest flood
 for mm=1:numel(stormHour)
  stormHourNearestFlood(mm)=hourPerDay*min(abs(stormHour(mm)-floodHour));
  % check if there's sea level data; if there are no coeval sea level data, set to nan
  if isnan(sl(find(dn==stormHour(mm))))  
   stormHourNearestFlood(mm)=nan;
  end
 end
 stormHourNearestFlood=stormHourNearestFlood';

 % reduce to daily values
 floodDay=unique(floor(floodHour));
 MslTideDay=unique(floor(tidalHour));
 floodDayWithStorm=nan*floodDay;
 tidalDay=zeros(size(floodDay));
 for mm=1:numel(floodDay)
  kk=[]; kk=find(floor(floodHour)==floodDay(mm));
  floodDayWithStorm(mm)=double(max(floodHourNearestStorm(kk)<=timeThreshold));
  tidalDay(mm)=double(max(floodTide(kk)));
 end

 stormDay=unique(floor(stormHour));
 stormDay(find(str2num(datestr(stormDay,10))>2016))=[];
 stormDayWithFlood=nan*floodDay;
 for mm=1:numel(stormDay)
  kk=[]; kk=find(floor(stormHour)==stormDay(mm));
  if ~isnan(sum(stormHourNearestFlood(kk)))
   stormDayWithFlood(mm)=double(max(stormHourNearestFlood(kk)<=timeThreshold));
  else
   stormDayWithFlood(mm)=nan;
  end
 end
 save(['noaa_tidegauge_',num2str(ID(nn)),'_ar_statistics_trsp',num2str(trspThreshold),'_dist',num2str(distThreshold),'.mat'],'floodDay*','stormDay*','tidalDay','MslTideDay','datum','sl','td','dn')
end

clearvars -except nnn
end
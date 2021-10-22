%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used to make Figures 2-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
load fileID.mat
distCrit=[1 2 1 2];
trspCrit=[500 500 250 250];

% load each gauge for each criteria
% compute: number of floods, number of storms, number of cooccurring floods and storms, fraction of floods that were storms
% for fraction, compute bootstrap errors that quantify uncertainty due to finite record length
% for fraction, also compute p value based on randomization give number of storms, floods, and data points

% define parameters and initialize some arrays
NID=numel(ID);
NCR=numel(distCrit);
NBT=1000; % 1,000 bootstrap iterations
NumFloods=zeros(NID,NCR,NBT);
NumStorms=zeros(NID,NCR,NBT);
NumFloSto=zeros(NID,NCR,NBT);
FracFloodsThatAreStorms=zeros(NID,NCR,NBT);
FracFloodsThatAreTides=zeros(NID,NCR,NBT);
FracFloodsThatAreStormsAndTides=zeros(NID,NCR,NBT);
FracStormsThatAreFloods=zeros(NID,NCR,NBT);
FracFloodsThatAreStormsRandom=zeros(NID,NCR,NBT);
rng(123);

for nn=1:NID, disp([num2str(nn)])
 for mm=1:NCR, disp(['...',num2str(mm)])
 %tic
  clearvars -except GT saveDir ID distCrit trspCrit nn mm Num* Frac* NID NCR NBT Lat Lon Name Threshold 
  load(['noaa_tidegauge_',num2str(ID(nn)),'_ar_statistics_trsp',num2str(trspCrit(mm)),'_dist',num2str(distCrit(mm)),'.mat'])
  floodDayWithTide=tidalDay; clear tidalDay
  Lat(nn)=datum.lat;
  Lon(nn)=datum.lon;
  Threshold(nn)=0.5+0.04*datum.GT; % sweet et al. 2018 threshold
  GT(nn)=datum.GT;
  Name(nn).name=datum.name;
  floodDayWithStormAndTide=floodDay(find(floodDayWithStorm==1&floodDayWithTide==1));
  floodDayWithStorm=floodDay(find(floodDayWithStorm==1));
  floodDayWithTide=floodDay(find(floodDayWithTide==1));
  stormDayWithFlood=stormDay(find(stormDayWithFlood==1));

  % create array of days during study period where you have 24 hourly water level values
  dn_days=[]; dn_days=datenum(1980,1,1):datenum(2016,12,31); 
  % these are the days between 1980-1-1 and 2016-12-31; only keep days with full data
  uu=[]; uu=sum(reshape(sl,24,numel(sl)/24),1);
  dn_days=dn_days(find(~isnan(uu)));
  clear uu
  nd=[]; nd=numel(dn_days);

  for kk=1:NBT
   % bootstrap data days
   ii=[]; ii=ceil(nd*rand(nd,1));

   NumMslTide(nn,mm,kk)=sum(ismember(dn_days(ii),MslTideDay));
   NumFloods(nn,mm,kk)=sum(ismember(dn_days(ii),floodDay));
   NumStorms(nn,mm,kk)=sum(ismember(dn_days(ii),stormDay));
   NumFloSto(nn,mm,kk)=sum(ismember(dn_days(ii),floodDayWithStorm));
   NumFloTid(nn,mm,kk)=sum(ismember(dn_days(ii),floodDayWithTide));
   NumFloStoTid(nn,mm,kk)=sum(ismember(dn_days(ii),floodDayWithStormAndTide));
 
   FracFloodsThatAreStorms(nn,mm,kk)=NumFloSto(nn,mm,kk)/NumFloods(nn,mm,kk);
   FracFloodsThatAreTides(nn,mm,kk)=NumFloTid(nn,mm,kk)/NumFloods(nn,mm,kk);
   FracFloodsThatAreStormsAndTides(nn,mm,kk)=NumFloStoTid(nn,mm,kk)/NumFloods(nn,mm,kk);
   FracStormsThatAreFloods(nn,mm,kk)=NumFloSto(nn,mm,kk)/NumStorms(nn,mm,kk);

   nm=[]; nm=NumMslTide(nn,mm,kk);
   nf=[]; nf=NumFloods(nn,mm,kk);
   ns=[]; ns=NumStorms(nn,mm,kk);

   % for clarity, simulate a daily poisson process
   % perform synthetic data significance test
   lambdastorm=[]; lambdaflood=[];
   lambdastorm=ns/nd; lambdaflood=nf/nd;
   x1=[]; x2=[];
   x1=poissrnd(lambdaflood,nd,1)';
   x2=poissrnd(lambdastorm,nd,1)';
   FracFloodsThatAreStormsRandom(nn,mm,kk)=(x1*x2')/sum(x1);
   FracStormsThatAreFloodsRandom(nn,mm,kk)=(x1*x2')/sum(x2);

  end
 %toc
 end
end

clearvars -except GT saveDir ID distCrit trspCrit nn mm Num* Frac* NID NCR NBT Lat Lon Name Threshold
[ll,ii]=sort(Lat);

ID=ID(ii)';
Threshold=Threshold(ii)';
Lat=Lat(ii)';
GT=GT(ii)';
Lon=Lon(ii)';
Name=Name(ii)';
FracFloodsThatAreTides=FracFloodsThatAreTides(ii,:,:);
FracFloodsThatAreStormsAndTides=FracFloodsThatAreStormsAndTides(ii,:,:);
FracFloodsThatAreStorms=FracFloodsThatAreStorms(ii,:,:);
FracFloodsThatAreStormsRandom=FracFloodsThatAreStormsRandom(ii,:,:);
FracStormsThatAreFloods=FracStormsThatAreFloods(ii,:,:);
FracStormsThatAreFloodsRandom=FracStormsThatAreFloodsRandom(ii,:,:);
NumMslTide=NumMslTide(ii,:,:);
NumStorms=NumStorms(ii,:,:);
NumFloods=NumFloods(ii,:,:);
NumFloSto=NumFloSto(ii,:,:);
NumFloTid=NumFloTid(ii,:,:);
NumFloStoTid=NumFloStoTid(ii,:,:);

clearvars -except GT ID Lat Lon Name Frac* Num* Threshold

% significance tests
sigTest1=mean(FracFloodsThatAreStorms>FracFloodsThatAreStormsRandom,3);
sigTest1=sigTest1>0.95;
sigTest1=double(sum(sigTest1,2)==4);

sigTest2=mean(FracStormsThatAreFloods>FracStormsThatAreFloodsRandom,3);
sigTest2=sigTest2>0.95;
sigTest2=double(sum(sigTest2,2)==4);

% start plotting
figure('color','white')

subplot(1,4,1)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,mean(NumFloods(:,1,:),3),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:50:200)
 caxis([0 200])
 title([{'a. Number'};{'HTF Days'}])
 set(gca,'fontweight','bold')

subplot(1,4,2)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,mean(NumStorms(:,1,:),3),'filled','markeredgecolor','k')
 set(gca,'fontweight','bold')
 colormap(turbo(12))
 c=colorbar('location','southoutside');
 set(c,'xtick',0:200:800)
 caxis([0 800])
 title([{'b. Number'};{'AR Days'}])
 set(gca,'fontweight','bold')

subplot(1,4,3)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,100*mean(FracFloodsThatAreStorms(:,1,:),3),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:10:40,'xticklabel',[{'0%'};{'10%'};{'20%'};{'30%'};{'40%'}])
 xtickangle(45)
 caxis([0 40])
 title([{'c. Percent HTF Days'};{'That Are AR Days'}])
 set(gca,'fontweight','bold')

 for nn=1:numel(Lat)
  if sigTest1(nn)==0
   plotm(Lat(nn),Lon(nn),'kx')
  end
 end

subplot(1,4,4)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,100*mean(FracStormsThatAreFloods(:,1,:),3),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:2.5:10,'xticklabel',[{'0%'};{'2.5%'};{'5%'};{'7.5%'};{'10%'}])
 xtickangle(45)
 caxis([0 10])
 title([{'d. Percent AR Days'};{'That Are HTF Days'}])
 set(gca,'fontweight','bold')

 for nn=1:numel(Lat)
  if sigTest2(nn)==0
   plotm(Lat(nn),Lon(nn),'kx')
  end
 end

chrisSaveFigurePng(['arsl_f02'],1000)

load clr0.mat
figure('color','white')
 subplot(5,1,1:2)
 hold on, box on, grid on
 for n=1:4
  plot(100*median(squeeze(FracFloodsThatAreStorms(:,n,:))'),'linewidth',2,'color',clr0(n,:))
 end
 for n=4:-1:1
  c=fill([1:24 24:-1:1],100*[prctile(squeeze(FracFloodsThatAreStorms(:,n,:))',97.5) fliplr(prctile(squeeze(FracFloodsThatAreStorms(:,n,:))',2.5))],[0.5 0.5 0.5]+0.5*clr0(n,:));
  set(c,'EdgeColor','none','facealpha',.25)
 end
  set(c,'EdgeColor','none','facealpha',.25)
 for n=4:-1:1
  plot(100*median(squeeze(FracFloodsThatAreStorms(:,n,:))'),'linewidth',2,'color',clr0(n,:))
 end
 legend([{'IVT \geq 500 kg/m/s; neareast 1 grid cell'};{'IVT \geq 500 kg/m/s; neareast 2 grid cells'};{'IVT \geq 250 kg/m/s; neareast 1 grid cell'};{'IVT \geq 250 kg/m/s; neareast 2 grid cells'}],'location','northeast','orientation','vertical')

 axis([1 24 0 100])
 set(gca,'xtick',1:24,'ytick',0:20:100,'yticklabel',[{'0%'};{'20%'};{'40%'};{'60%'};{'80%'};{'100%'}])
 title([{'a. Percent HTF Days That Are AR Days'}])
 set(gca,'fontweight','bold')
 set(gca,'xticklabel',[{'San Diego'};{'La Jolla'};{'Los Angeles'};{'Santa Monica'};{'Port San Luis'};{'Monterey'};{'Alameda'};{'San Francisco'};{'Point Reyes'};{'Port Chicago'};{'Arena Cove'};{'Humboldt Bay'};{'Crescent City'};{'Port Orford'};{'Charleston'};{'South Beach'};{'Astoria'};{'Toke Point'};{'Seattle'};{'Port Townsend'};{'Port Angeles'};{'Neah Bay'};{'Friday Harbor'};{'Cherry Point'}])

clear all

load fileID.mat
distCrit=1;
trspCrit=500;
NID=numel(ID);
NCR=numel(distCrit);
NBT=1000; % 1,000 bootstrap iterations

% first create structures
rng(123);
%dn_days=datenum(1980,1,1):datenum(2016,12,31);
%yr_days=str2num(datestr(dn_days,10));
% use meteorological years as in sweet et al. noaa htf reports
dn_days=datenum(1980,5,1):datenum(2016,4,30);
yr_days=str2num(datestr(dn_days,10));
mo_days=str2num(datestr(dn_days,5));
% if month is 1-4 make it previous year
yr_days(find(mo_days<=4))=yr_days(find(mo_days<=4))-1;

NumMslTide=nan*ones(NID,numel(dn_days));
NumFloods=nan*ones(NID,numel(dn_days));
NumStorms=nan*ones(NID,numel(dn_days));
NumTides=nan*ones(NID,numel(dn_days));

for nn=1:NID, disp([num2str(nn)])
 clearvars -except GT saveDir ID distCrit trspCrit nn mm Num* Frac* NID NCR NBT Lat Lon Name Threshold *_days
 load(['noaa_tidegauge_',num2str(ID(nn)),'_ar_statistics_trsp',num2str(trspCrit),'_dist',num2str(distCrit),'.mat'])
 NumMslTide(nn,:)=ismember(dn_days,MslTideDay);
 NumTides(nn,:)=ismember(dn_days,floodDay(find(tidalDay==1)));
 NumFloods(nn,:)=ismember(dn_days,floodDay);
 NumStorms(nn,:)=ismember(dn_days,stormDay);
end

clearvars -except NID NCR NBT *_days Num*
nd=numel(dn_days);
ns=NID;
years=1980:2015;

MslTideYear=nan*zeros(numel(years),NBT);
FloodYear=nan*zeros(numel(years),NBT);
StormYear=nan*zeros(numel(years),NBT);
TideYear=nan*zeros(numel(years),NBT);
Msl=nan*zeros(numel(years),NBT);
FloodYearPoisson=nan*zeros(numel(years),NBT);
StormYearPoisson=nan*zeros(numel(years),NBT);
MslRandom=nan*zeros(numel(years),NBT);

load('/Users/christopherpiecuch/Documents/Data/atmospheric_rivers/20210817_sealevel_statistics_w_precip_10.mat')
arraySL=nan*ones(NID,numel(dn_days));
for kk=1:NID
 i1=[];i2=[];i1=find(data(kk).date==min(dn_days));i2=find(data(kk).date==max(dn_days));
 arraySL(kk,:)=data(kk).sl_dt_ds(i1:i2)-data(kk).sl_dt_ds_hf(i1:i2);
end
clear data i1 i2

for kk=1:NBT, disp(num2str(kk))
 ii=[]; ii=ceil(ns*rand(ns,1));
 jj=[]; jj=ceil(nd*rand(nd,1));
 NumMslTideTemp=[]; NumMslTideTemp=NumMslTide(ii,jj);
 NumTidesTemp=[]; NumTidesTemp=NumTides(ii,jj);
 NumFloodsTemp=[]; NumFloodsTemp=NumFloods(ii,jj);
 NumStormsTemp=[]; NumStormsTemp=NumStorms(ii,jj);
 arraySLTemp=[]; arraySLTemp=arraySL(ii,jj);
 yr_daysTemp=[]; yr_daysTemp=yr_days(jj);
 for ll=1:numel(years)
  uu=[]; uu=find(yr_daysTemp==years(ll));
  TideYear(ll,kk)=sum(sum(NumTidesTemp(:,uu)))/ns;
  FloodYear(ll,kk)=sum(sum(NumFloodsTemp(:,uu)))/ns;
  MslTideYear(ll,kk)=sum(sum(NumMslTideTemp(:,uu)))/ns;
  StormYear(ll,kk)=sum(sum(NumStormsTemp(:,uu)))/ns;
  Msl(ll,kk)=nanmean(nanmean(arraySLTemp(:,uu)));
 end
 % now do stochastic simulation
 lambdaFlood=[]; lambdaStorm=[];
 lambdaFlood=sum(FloodYear(:,kk))/numel(years);
 lambdaStorm=sum(StormYear(:,kk))/numel(years);
 FloodYearPoisson(:,kk)=poissrnd(lambdaFlood,numel(years),1);
 StormYearPoisson(:,kk)=poissrnd(lambdaStorm,numel(years),1);
 MslRandom(:,kk)=randn(size(Msl(:,kk)))*std(Msl(:,kk));

 CorrFloodTide(kk)=corr(FloodYear(:,kk),TideYear(:,kk),'Type','Pearson');
 CorrFloodStorm(kk)=corr(FloodYear(:,kk),StormYear(:,kk),'Type','Pearson');
 CorrFloodStormPoisson(kk)=corr(FloodYearPoisson(:,kk),StormYearPoisson(:,kk),'Type','Pearson');
 CorrFloodMsl(kk)=corr(FloodYear(:,kk),Msl(:,kk),'Type','Pearson');
 CorrFloodMslTide(kk)=corr(FloodYear(:,kk),MslTideYear(:,kk),'Type','Pearson');
 CorrFloodMslRandom(kk)=corr(FloodYearPoisson(:,kk),MslRandom(:,kk),'Type','Pearson');
 AutoCorrMsl(kk)=autocorrelation(Msl(:,kk),1);

end

% here is where you shade things in
load clr0.mat
ax1=subplot(5,1,4:5);

yyaxis left
plot(years,prctile(FloodYear',50),'-','linewidth',2,'color',clr0(1,:))
hold on
plot(years,prctile(FloodYear',50),'--','linewidth',2,'color',clr0(1,:))
plot(years,prctile(FloodYear',50),'-','linewidth',2,'color',clr0(2,:))
plot(years,prctile(FloodYear',50),'-','linewidth',2,'color','k')
plot(years,prctile(MslTideYear',50),'--','linewidth',2,'color',clr0(1,:))
c=fill([years fliplr(years)],[prctile(FloodYear',97.5) fliplr(prctile(FloodYear',2.5))],[0.5 0.5 0.5]+0.5*clr0(1,:));
set(c,'EdgeColor','none','facealpha',.25)
hold on
plot(years,prctile(FloodYear',50),'-','linewidth',2,'color',clr0(1,:))
plot(years,prctile(MslTideYear',50),'--','linewidth',2,'color',clr0(1,:))
pause(0.1)
grid on, box on, hold on
set(gca,'xtick',1980:5:2015,'ytick',0:5:20,'fontweight','bold')
ylim([0 20]) 
pause(0.1)
ylabel('HTF Days')
xlabel('Meteorological Year')
title('b. Annual HTF Days, AR Days, and Mean Sea Level','fontweight','bold')

yyaxis right
c=fill([years fliplr(years)],[prctile(StormYear',97.5) fliplr(prctile(StormYear',2.5))],[0.5 0.5 0.5]+0.5*clr0(2,:));
set(c,'EdgeColor','none','facealpha',.25)
hold on
plot(years,prctile(StormYear',50),'-','linewidth',2,'color',clr0(2,:))
ylim([0 40])
pause(0.1)
grid on, box on, hold on
set(gca,'xtick',1980:5:2015,'ytick',0:10:40,'fontweight','bold')
ylabel('AR Days')
legend([{'Observed HTF days'};{'Hypothetical HTF days'};{'AR days'};{'MSL'}],'location','northwest','orientation','horizontal')

ax2 = axes('position', ax1.Position);
c=fill([years fliplr(years)],1e2*[prctile(Msl',97.5) fliplr(prctile(Msl',2.5))],[0.5 0.5 0.5]);
set(c,'EdgeColor','none','facealpha',.25)
hold on
plot(years,1e2*prctile(Msl',50),'-','linewidth',2,'color','k')
ax2.Color='none'; 
grid(ax2,'on')
grid on, box on, hold on
set(gca,'xtick',1980:5:2015,'ytick',-5:5:15,'fontweight','bold')
ylim([-5 15])
ax2.YTickLabel = strcat(ax2.YTickLabel, {'       '}, {'       '}); 
ylabel('Mean Sea Level (cm)')

chrisSaveFigurePng(['arsl_f03'],1000)
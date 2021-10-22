%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used to make Figures 1, 4, 5, and S1-S3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define parameters, initialize arrays
% and load composite analysis data
load('ar_and_ss_statistics_w_precip_10.mat')
rng(123);
NID=numel(arss);
NBT=1000; % 1,000 bootstrap iterations
lambda=0.00:0.01:1; % values of ridge parameter

% perform regression analysis (equation 1 in the main text)
% solve regression model three ways:
% least squares (LS), ridge regression (RR),
% and james stein estimation (JS)
RegCoefLS=nan*zeros(8,NID,NBT);
RegCoefJS=nan*zeros(8,NID,NBT);
RegCoefRR=nan*zeros(8,numel(lambda),NID,NBT);
VarExplModelLS=nan*zeros(NID,NBT);
VarExplModelJS=nan*zeros(NID,NBT);
VarExplModelRR=nan*zeros(numel(lambda),NID,NBT);
VarExplPresLS=nan*zeros(NID,NBT);
VarExplTauxLS=nan*zeros(NID,NBT);
VarExplTauyLS=nan*zeros(NID,NBT);
VarExplPrecipLS=nan*zeros(NID,NBT);
VarExplPresJS=nan*zeros(NID,NBT);
VarExplTauxJS=nan*zeros(NID,NBT);
VarExplTauyJS=nan*zeros(NID,NBT);
VarExplPrecipJS=nan*zeros(NID,NBT);
VarExplPresRR=nan*zeros(numel(lambda),NID,NBT);
VarExplTauxRR=nan*zeros(numel(lambda),NID,NBT);
VarExplTauyRR=nan*zeros(numel(lambda),NID,NBT);
VarExplPrecipRR=nan*zeros(numel(lambda),NID,NBT);
MeanModelRR=zeros(NID,NBT,numel(lambda));
StdModelRR=zeros(NID,NBT,numel(lambda));
MeanPresRR=zeros(NID,NBT,numel(lambda));
MeanTauxRR=zeros(NID,NBT,numel(lambda));
MeanTauyRR=zeros(NID,NBT,numel(lambda));
MeanPrecipRR=zeros(NID,NBT,numel(lambda));

for nn=1:NID, disp(num2str(nn))

 % define longitude and latitude as well as number of events
 Lat(nn)=arss(nn).latitude;
 Lon(nn)=arss(nn).longitude;
 NN=[]; NN=numel(arss(nn).data);

 % perform analysis using bootstrapping
 for mm=1:NBT
  ii=[]; ii=ceil(NN*rand(NN,1)); % bootstrap indices
  % fit model
  y=[]; X=[]; S=[]; betaBS=[]; betaLS=[]; betaJS=[]; betaRR=[];

  % define vector of observations y and design matrix X
  y=arss(nn).data(ii);
  X=[arss(nn).pres(ii) arss(nn).Hpres(ii) arss(nn).taux(ii) arss(nn).Htaux(ii) arss(nn).tauy(ii) arss(nn).Htauy(ii) arss(nn).precip(ii) arss(nn).Hprecip(ii)];

  % baseline determined using Matlab's regression call
  betaBS=regress(y,[X ones(NN,1)]);
  betaBS(numel(betaBS))=[]; % delete regression on constant

  % do transformations needed for james-stein estimation and ridge regression
  ybar=mean(y);
  y=y-ybar;
  Xbar=mean(X);
  X=X-ones(NN,1)*Xbar;
  Xstd=sqrt(NN-1)*std(X);
  X=X./(ones(NN,1)*Xstd);
  S=X'*X;

  % solve regression problem using least squares, james-stein, and ridge regression
  betaLS=inv(S)*X'*y;
  betaJS=(1-(((numel(X)/NN)-2)*var(y-X*betaLS))/(betaLS'*S*betaLS))*betaLS;
  for kk=1:numel(lambda), %disp(num2str(kk))
   betaRR(:,kk)=inv(S+lambda(kk)*eye(size(S)))*X'*y;
  end

  % scale back to units
  betaLS=betaLS./Xstd'; % scale back to units
  betaJS=betaJS./Xstd'; % scale back to units
  betaRR=betaRR./(Xstd'*ones(1,numel(lambda))); % scale back to units
  % store
  RegCoefLS(:,nn,mm)=betaLS;
  RegCoefJS(:,nn,mm)=betaJS;
  RegCoefRR(:,:,nn,mm)=betaRR;

  % redefine vector of observations y and design matrix X
  y=arss(nn).data(ii);
  X=[arss(nn).pres(ii) arss(nn).Hpres(ii) arss(nn).taux(ii) arss(nn).Htaux(ii) arss(nn).tauy(ii) arss(nn).Htauy(ii) arss(nn).precip(ii) arss(nn).Hprecip(ii)];

  MeanData(nn,mm)=mean(y);
  StdData(nn,mm)=std(y);
  for kk=1:numel(lambda), %disp(num2str(kk))
   MeanModelRR(nn,mm,kk)=mean(X*betaRR(:,kk));
   StdModelRR(nn,mm,kk)=std(X*betaRR(:,kk));
   MeanPresRR(nn,mm,kk)=mean(X(:,1:2)*betaRR(1:2,kk));
   MeanTauxRR(nn,mm,kk)=mean(X(:,3:4)*betaRR(3:4,kk));
   MeanTauyRR(nn,mm,kk)=mean(X(:,5:6)*betaRR(5:6,kk));
   MeanPrecipRR(nn,mm,kk)=mean(X(:,7:8)*betaRR(7:8,kk));
  end

  MeanModel(nn,mm)=mean(X*betaLS);
  StdModel(nn,mm)=std(X*betaLS);
  MeanPres(nn,mm)=mean(X(:,1:2)*betaLS(1:2));
  MeanTaux(nn,mm)=mean(X(:,3:4)*betaLS(3:4));
  MeanTauy(nn,mm)=mean(X(:,5:6)*betaLS(5:6));
  MeanPrecip(nn,mm)=mean(X(:,7:8)*betaLS(7:8));

  VarExplModelLS(nn,mm)=1-var(y-X*betaLS)/var(y);
  VarExplModelJS(nn,mm)=1-var(y-X*betaJS)/var(y);
  for kk=1:numel(lambda), %disp(num2str(kk))
   VarExplModelRR(kk,nn,mm)=1-var(y-X*betaRR(:,kk))/var(y);
  end

  VarExplPresJS(nn,mm)=1-var(y-X(:,1:2)*betaJS(1:2))/var(y);
  VarExplTauxJS(nn,mm)=1-var(y-X(:,3:4)*betaJS(3:4))/var(y);
  VarExplTauyJS(nn,mm)=1-var(y-X(:,5:6)*betaJS(5:6))/var(y);
  VarExplPrecipJS(nn,mm)=1-var(y-X(:,7:8)*betaJS(7:8))/var(y);

  VarExplPresLS(nn,mm)=1-var(y-X(:,1:2)*betaLS(1:2))/var(y);
  VarExplTauxLS(nn,mm)=1-var(y-X(:,3:4)*betaLS(3:4))/var(y);
  VarExplTauyLS(nn,mm)=1-var(y-X(:,5:6)*betaLS(5:6))/var(y);
  VarExplPrecipLS(nn,mm)=1-var(y-X(:,7:8)*betaLS(7:8))/var(y);

  for kk=1:numel(lambda), %disp(num2str(kk))
   VarExplPresRR(kk,nn,mm)=1-var(y-X(:,1:2)*betaRR(1:2,kk))/var(y);
   VarExplTauxRR(kk,nn,mm)=1-var(y-X(:,3:4)*betaRR(3:4,kk))/var(y);
   VarExplTauyRR(kk,nn,mm)=1-var(y-X(:,5:6)*betaRR(5:6,kk))/var(y);
   VarExplPrecipRR(kk,nn,mm)=1-var(y-X(:,7:8)*betaRR(7:8,kk))/var(y);
  end

 end
end

% perform bootstrapping with constant coefficients
lambdaFixed=0.3; % ridge regression lambda value

for mm=1:NBT
 clear y X S I betaRR ii y x *Store ybar Xstd i2u
 y=[]; X=[]; S=[]; I=[]; betaRR=[];
 for nn=1:NID
  NN=[]; NN=numel(arss(nn).data);
  ii=[]; ii=ceil(NN*rand(NN,1)); % bootstrap indices
  I((numel(I)+1):(numel(I)+numel(ii)))=nn;
  y=[y; arss(nn).data(ii)];
  x=[]; x=[arss(nn).pres(ii) arss(nn).Hpres(ii) arss(nn).taux(ii) arss(nn).Htaux(ii) arss(nn).tauy(ii) arss(nn).Htauy(ii) arss(nn).precip(ii) arss(nn).Hprecip(ii)];
  X=[X; x];
 end

  % do transformations needed for james-stein estimation and ridge regression
  yStore=[]; yStore=y;
  XStore=[]; XStore=X;
  ybar=mean(y);
  y=y-ybar;
  Xbar=mean(X);
  X=X-ones(numel(I),1)*Xbar;
  Xstd=sqrt(numel(I)-1)*std(X);
  X=X./(ones(numel(I),1)*Xstd);
  S=X'*X;

  % solve regression problem using least squares, james-stein, and ridge regression
  betaRR=inv(S+lambdaFixed*eye(size(S)))*X'*y;

  % scale back to units
  betaRR=betaRR./Xstd'; % scale back to units
  
  % return values
  y=yStore; X=XStore;
 
  for nn=1:NID
   i2u=[]; i2u=find(I==nn);
   VarExplModelRRConst(nn,mm)=1-var(y(i2u)-X(i2u,:)*betaRR)/var(y(i2u));
   VarExplPresRRConst(nn,mm)=1-var(y(i2u)-X(i2u,1:2)*betaRR(1:2))/var(y(i2u));
   VarExplTauxRRConst(nn,mm)=1-var(y(i2u)-X(i2u,3:4)*betaRR(3:4))/var(y(i2u));
   VarExplTauyRRConst(nn,mm)=1-var(y(i2u)-X(i2u,5:6)*betaRR(5:6))/var(y(i2u));
   VarExplPrecipRRConst(nn,mm)=1-var(y(i2u)-X(i2u,7:8)*betaRR(7:8))/var(y(i2u));
   MeanModelRRConst(nn,mm)=mean(X(i2u,:)*betaRR);
   StdModelRRConst(nn,mm)=std(X(i2u,:)*betaRR);
   MeanTauxRRConst(nn,mm)=mean(X(i2u,3:4)*betaRR(3:4));
   MeanTauyRRConst(nn,mm)=mean(X(i2u,5:6)*betaRR(5:6));
   MeanPresRRConst(nn,mm)=mean(X(i2u,1:2)*betaRR(1:2));
   MeanPrecipRRConst(nn,mm)=mean(X(i2u,7:8)*betaRR(7:8));
  end
end

for nn=1:NID
 Name(nn).name=arss(nn).name;
end

clearvars -except NBT Atm* Mean* Std* Var* Reg* Lon Lat arss Name

% sort all data by latitude for plotting purposes
[ll,ii]=sort(Lat);
Name=Name(ii);
Lat=Lat(ii)';
Lon=Lon(ii)';
RegCoefJS=RegCoefJS(:,ii,:);
RegCoefLS=RegCoefLS(:,ii,:);
RegCoefRR=RegCoefRR(:,:,ii,:);
RegCoefPrecipRRa=squeeze(RegCoefRR(7,:,:,:));
RegCoefPrecipRRb=squeeze(RegCoefRR(8,:,:,:));
RegCoefPresRRa=squeeze(RegCoefRR(1,:,:,:));
RegCoefPresRRb=squeeze(RegCoefRR(2,:,:,:));
RegCoefTauxRRa=squeeze(RegCoefRR(3,:,:,:));
RegCoefTauxRRb=squeeze(RegCoefRR(4,:,:,:));
RegCoefTauyRRa=squeeze(RegCoefRR(5,:,:,:));
RegCoefTauyRRb=squeeze(RegCoefRR(6,:,:,:));

MeanData=MeanData(ii,:);
StdData=StdData(ii,:);

MeanModel=MeanModel(ii,:);
MeanPrecip=MeanPrecip(ii,:);
MeanPres=MeanPres(ii,:);
MeanTaux=MeanTaux(ii,:);
MeanTauy=MeanTauy(ii,:);

MeanModelRR=MeanModelRR(ii,:,:);
MeanPrecipRR=MeanPrecipRR(ii,:,:);
MeanPresRR=MeanPresRR(ii,:,:);
MeanTauxRR=MeanTauxRR(ii,:,:);
MeanTauyRR=MeanTauyRR(ii,:,:);

StdModel=StdModel(ii,:);
StdModelRR=StdModelRR(ii,:,:);

MeanModelRRConst=MeanModelRRConst(ii,:);
MeanPrecipRRConst=MeanPrecipRRConst(ii,:);
MeanPresRRConst=MeanPresRRConst(ii,:);
MeanTauxRRConst=MeanTauxRRConst(ii,:);
MeanTauyRRConst=MeanTauyRRConst(ii,:);
StdModelRRConst=StdModelRRConst(ii,:);

VarExplModelLS=VarExplModelLS(ii,:);
VarExplTauyLS=VarExplTauyLS(ii,:);
VarExplPresLS=VarExplPresLS(ii,:);
VarExplPrecipLS=VarExplPrecipLS(ii,:);
VarExplTauxLS=VarExplTauxLS(ii,:);

VarExplModelJS=VarExplModelJS(ii,:);
VarExplTauyJS=VarExplTauyJS(ii,:);
VarExplPresJS=VarExplPresJS(ii,:);
VarExplPrecipJS=VarExplPrecipJS(ii,:);
VarExplTauxJS=VarExplTauxJS(ii,:);

VarExplModelRR=permute(VarExplModelRR,[2 3 1]);
VarExplTauyRR=permute(VarExplTauyRR,[2 3 1]);
VarExplPresRR=permute(VarExplPresRR,[2 3 1]);
VarExplPrecipRR=permute(VarExplPrecipRR,[2 3 1]);
VarExplTauxRR=permute(VarExplTauxRR,[2 3 1]);
VarExplModelRR=VarExplModelRR(ii,:,:);
VarExplTauyRR=VarExplTauyRR(ii,:,:);
VarExplPresRR=VarExplPresRR(ii,:,:);
VarExplPrecipRR=VarExplPrecipRR(ii,:,:);
VarExplTauxRR=VarExplTauxRR(ii,:,:);

VarExplModelRRConst=VarExplModelRRConst(ii,:);
VarExplTauyRRConst=VarExplTauyRRConst(ii,:);
VarExplPresRRConst=VarExplPresRRConst(ii,:);
VarExplPrecipRRConst=VarExplPrecipRRConst(ii,:);
VarExplTauxRRConst=VarExplTauxRRConst(ii,:);
clearvars -except NBT Atm* Mean* Std* Var* Reg* Lon Lat  arss Name

% start making figures

lamToUse=31; % lambda=0.3; 31st index of 101-element vector 0:0.01:1

% Figure 4
figure('color','white')

subplot(1,4,1)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,1e2*mean(MeanData'),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:5:25)
 caxis([0 20])
 title([{'a. Observed Mean'};{'Storm Surge (cm)'};])
 set(gca,'fontweight','bold')

subplot(1,4,2)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,1e2*mean(StdData'),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:5:25)
 caxis([0 20])
 title([{'b. Observed Surge'};{'Standard Deviation (cm)'};])
 set(gca,'fontweight','bold')

subplot(1,4,3)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,1e2*mean(squeeze(MeanModelRR(:,:,lamToUse))'),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:5:25)
 caxis([0 20])
 title([{'c. Modeled Mean'};{'Storm Surge (cm)'};])
 set(gca,'fontweight','bold')

subplot(1,4,4)
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 setm(gca,'parallellabel','off','meridianlabel','off')
 scatterm(Lat,Lon,50,1e2*mean(squeeze(StdModelRR(:,:,lamToUse))'),'filled','markeredgecolor','k')
 colormap(turbo(12))
 set(gca,'fontweight','bold')
 c=colorbar('location','southoutside');
 set(c,'xtick',0:5:25)
 caxis([0 20])
 title([{'d. Modeled Surge'};{'Standard Deviation (cm)'};])
 set(gca,'fontweight','bold')

chrisSaveFigurePng(['arsl_f04'],1000)

load clr0.mat
close all
siteNames=[{'San Diego'};{'La Jolla'};{'Los Angeles'};{'Santa Monica'};{'Port San Luis'};{'Monterey'};{'Alameda'};{'San Francisco'};{'Point Reyes'};{'Port Chicago'};{'Arena Cove'};{'Humboldt Bay'};{'Crescent City'};{'Port Orford'};{'Charleston'};{'South Beach'};{'Astoria'};{'Toke Point'};{'Seattle'};{'Port Townsend'};{'Port Angeles'};{'Neah Bay'};{'Friday Harbor'};{'Cherry Point'}];
letters='abcd';
clr0(1,:)=0*clr0(1,:);

figure('color','white')
 subplot(5,1,1:2)
 hold on, box on, grid on
 plot(1e2*median(squeeze(MeanData)'),'linewidth',2,'color',clr0(1,:))
 plot(1e2*median(squeeze(MeanModelRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(2,:))
 plot(1e2*median(squeeze(MeanTauxRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(3,:))
 plot(1e2*median(squeeze(MeanTauyRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(4,:))
 plot(1e2*median(squeeze(MeanPresRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(5,:))
 plot(1e2*median(squeeze(MeanPrecipRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(6,:))

  c=fill([1:24 24:-1:1],1e2*[prctile(MeanData',97.5) fliplr(prctile(MeanData',2.5))],[0.5 0.5 0.5]+0.5*clr0(1,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],1e2*[prctile(squeeze(MeanModelRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(MeanModelRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(2,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],1e2*[prctile(squeeze(MeanTauxRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(MeanTauxRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(3,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],1e2*[prctile(squeeze(MeanTauyRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(MeanTauyRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(4,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],1e2*[prctile(squeeze(MeanPresRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(MeanPresRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(5,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],1e2*[prctile(squeeze(MeanPrecipRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(MeanPrecipRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(6,:)); set(c,'EdgeColor','none','facealpha',.25)

 plot(1e2*median(squeeze(MeanData)'),'linewidth',2,'color',clr0(1,:))
 plot(1e2*median(squeeze(MeanModelRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(2,:))
 plot(1e2*median(squeeze(MeanTauxRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(3,:))
 plot(1e2*median(squeeze(MeanTauyRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(4,:))
 plot(1e2*median(squeeze(MeanPresRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(5,:))
 plot(1e2*median(squeeze(MeanPrecipRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(6,:))

 axis([1 24 -5 20])
 set(gca,'xtick',1:24,'ytick',-5:5:25,'xticklabel',[])
 title('a. Observed and Modeled Mean Storm Surge')
 legend boxoff
 set(gca,'fontweight','bold')
 ylabel('Surge (cm)','fontweight','bold')

 legend([{'$\zeta$'};{'$\hat{\zeta}$'};{'$\hat{\zeta}_{\pi}$'};{'$\hat{\zeta}_{\tau}$'};{'$\hat{\zeta}_{p}$'};{'$\hat{\zeta}_{q}$'}],'location','northwest','orientation','horizontal','Interpreter','latex')

 subplot(5,1,3:4)

 hold on, box on, grid on

 plot(100*median(squeeze(VarExplModelRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(2,:))
 plot(100*median(squeeze(VarExplTauxRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(3,:))
 plot(100*median(squeeze(VarExplTauyRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(4,:))
 plot(100*median(squeeze(VarExplPresRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(5,:))
 plot(100*median(squeeze(VarExplPrecipRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(6,:))

  c=fill([1:24 24:-1:1],100*[prctile(squeeze(VarExplModelRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(VarExplModelRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(2,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],100*[prctile(squeeze(VarExplTauxRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(VarExplTauxRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(3,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],100*[prctile(squeeze(VarExplTauyRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(VarExplTauyRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(4,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],100*[prctile(squeeze(VarExplPresRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(VarExplPresRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(5,:)); set(c,'EdgeColor','none','facealpha',.25)
  c=fill([1:24 24:-1:1],100*[prctile(squeeze(VarExplPrecipRR(:,:,lamToUse))',97.5) fliplr(prctile(squeeze(VarExplPrecipRR(:,:,lamToUse))',2.5))],[0.5 0.5 0.5]+0.5*clr0(6,:)); set(c,'EdgeColor','none','facealpha',.25)

 plot(100*median(squeeze(VarExplModelRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(2,:))
 plot(100*median(squeeze(VarExplTauxRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(3,:))
 plot(100*median(squeeze(VarExplTauyRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(4,:))
 plot(100*median(squeeze(VarExplPresRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(5,:))
 plot(100*median(squeeze(VarExplPrecipRR(:,:,lamToUse))'),'linewidth',2,'color',clr0(6,:))

 axis([1 24 0 100])
 set(gca,'xtick',1:24,'ytick',0:20:100,'yticklabel',[{'0%'};{'20%'};{'40%'};{'60%'};{'80%'};{'100%'}])
 title('b. Percent of Observed Storm-Surge Variance Explained by Model')
 set(gca,'fontweight','bold')
 ylabel('Percent','fontweight','bold')

 set(gca,'xticklabel',[{'San Diego'};{'La Jolla'};{'Los Angeles'};{'Santa Monica'};{'Port San Luis'};{'Monterey'};{'Alameda'};{'San Francisco'};{'Point Reyes'};{'Port Chicago'};{'Arena Cove'};{'Humboldt Bay'};{'Crescent City'};{'Port Orford'};{'Charleston'};{'South Beach'};{'Astoria'};{'Toke Point'};{'Seattle'};{'Port Townsend'};{'Port Angeles'};{'Neah Bay'};{'Friday Harbor'};{'Cherry Point'}])

chrisSaveFigurePng(['arsl_f05'],1000)

load clr0.mat

load('theoretical_coefficients','a*','b*')
figure('color','white')
LAM=0:0.01:1; NLAM=numel(LAM);

subplot(4,2,1)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(ax); yy=max(ax);
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefTauxRRa;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1,'xticklabel',[])
title('$a_{\pi}$','Interpreter','latex','fontsize',16)
axis([0 1 -0.4 0.4])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('m/Pa')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

subplot(4,2,2)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(bx); yy=max(bx);
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefTauxRRb;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1,'xticklabel',[])
title('$b_{\pi}$','Interpreter','latex','fontsize',16)
axis([0 1 -0.4 0.4])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('m/Pa')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

subplot(4,2,3)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(ay); yy=max(ay);
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefTauyRRa;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1,'xticklabel',[])
title('$a_{\tau}$','Interpreter','latex','fontsize',16)
legend boxoff
axis([0 1 -0.4 0.4])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('m/Pa')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')
legend([{'Theoretical'};{'Empirical'}],'location','southeast','orientation','horizontal')

subplot(4,2,4)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(by); yy=max(by);
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefTauyRRb;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1,'xticklabel',[])
title('$b_{\tau}$','Interpreter','latex','fontsize',16)
axis([0 1 -0.4 0.4])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('m/Pa')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

subplot(4,2,5)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(ap); yy=max(ap);
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefPresRRa;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1,'xticklabel',[])
title('$a_{p}$','Interpreter','latex','fontsize',16)
axis([0 1 -2e-4 0])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('m/Pa')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

subplot(4,2,6)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(bp); yy=max(bp);
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefPresRRb;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1,'xticklabel',[])
title('$b_{p}$','Interpreter','latex','fontsize',16)
axis([0 1 -1e-4 1e-4])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('m/Pa')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

subplot(4,2,7)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(aq)/86400; yy=max(aq)/86400;
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefPrecipRRa/86400;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1)
title('$a_{q}$','Interpreter','latex','fontsize',16)
set(gca,'xtick',0:0.2:1)
xlabel('$\lambda$','Interpreter','latex','fontsize',16)
axis([0 1 -2 6])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('days')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

subplot(4,2,8)
hold on, box on, grid on
xx=[]; yy=[];
xx=min(bq)/86400; yy=max(bq)/86400;
c=fill([0 1 1 0 0],[xx xx yy yy xx],clr0(1,:));
set(c,'edgecolor',clr0(1,:),'facecolor',clr0(1,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
XX=RegCoefPrecipRRb/86400;
XX=reshape(XX,NLAM,numel(Lat)*NBT);
xx=[]; yy=[];
xx=prctile(XX',2.5); yy=fliplr(prctile(XX',97.5));
c=fill([LAM fliplr(LAM)],[xx yy],clr0(2,:));
set(c,'edgecolor',clr0(2,:),'facecolor',clr0(2,:),'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
set(gca,'xtick',0:0.2:1)
title('$b_{q}$','Interpreter','latex','fontsize',16)
set(gca,'xtick',0:0.2:1)
xlabel('$\lambda$','Interpreter','latex','fontsize',16)
axis([0 1 -2 6])
AX=[]; AX=axis; set(gca,'ytick',AX(3):((AX(4)-AX(3))/4):AX(4))
grid on, box on, hold on
ylabel('days')
set(gca,'fontsize',12)
plot([0.3 0.3],[-10 10],'k--')

chrisSaveFigurePng(['arsl_fS2'],1000)
close all

figure('color','white')
LMNO=2:2:24;
for nn=1:9
 subplot(3,3,nn)
 histogram(arss(LMNO(nn)).data*100,[-10:5:100],'Normalization','probability','edgecolor',clr0(1,:),'facecolor',[.5 .5 .5]+0.5*clr0(1,:))
 axis([0 100 0 0.4])
 set(gca,'fontsize',10,'fontweight','bold','xtick',0:20:100,'ytick',0:0.1:0.4)
 hold on, box on, grid on
 c=fill([56 64 64 56 56],[0 0 1 1 0],[.5 .5 .5]);
 set(c,'facealpha',0.25,'edgealpha',0.25,'linewidth',2)
 title([arss(LMNO(nn)).name,' (n=',num2str(numel(arss(LMNO(nn)).data)),')'],'fontsize',10)
 xlabel('Surge (cm)','fontsize',10)
 ylabel('Probability','fontsize',10)
end

chrisSaveFigurePng(['arsl_fS1'],1000)
close all

figure('color','white')
LMNO=2:2:24;
for nn=1:9
 subplot(3,3,nn)
 plot([0 2],[0 2],'--','linewidth',2,'color',[.5 .5 .5])
 hold on
 plot(abs(arss(LMNO(nn)).taux),abs(arss(LMNO(nn)).tauy),'o','markerfacecolor',clr0(1,:),'markeredgecolor','k','markersize',5)
 zzz=[]; zzz=((mean(abs(arss(LMNO(nn)).taux)<abs(arss(LMNO(nn)).tauy))));
 axis([0 2 0 2])
 set(gca,'fontsize',10,'fontweight','bold','xtick',0:0.5:2,'ytick',0:0.5:2)
 hold on, box on, grid on
 title([arss(LMNO(nn)).name,' (',num2str(round(1e2*zzz)),'%)'],'fontsize',10)
 xlabel('| \pi | (N/m^2)','fontsize',10)
 ylabel('| \tau | (N/m^2)','fontsize',10)
end

chrisSaveFigurePng(['arsl_fS3'],1000)
close all

clearvars -except NBT Lon Lat 
load('/Users/christopherpiecuch/Documents/Data/atmospheric_rivers/SIO_R1_1948-2017_Comprehensive.mat','Longitude','Latitude')
uniqueCoords=unique(Longitude+i*Latitude);
clear Longitude Latitude
Longitude=real(uniqueCoords);
Latitude=imag(uniqueCoords);
clear uniqueCoords

figure('color','white')
 worldmap([min(Lat)-1 max(Lat)+1],[min(Lon)-1 max(Lon)+1])
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.7 .7 .7])
 geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 for n=1:numel(Longitude)
  plotm(Latitude(n),Longitude(n),'s','markersize',5,'markeredgecolor','k','markerfacecolor','k')
  plotm((Latitude(n)+1.25)*ones(1,26),(Longitude(n)-1.25):0.1:(Longitude(n)+1.25),'k')
  plotm((Latitude(n)-1.25)*ones(1,26),(Longitude(n)-1.25):0.1:(Longitude(n)+1.25),'k')
  plotm((Latitude(n)-1.25):0.1:(Latitude(n)+1.25),(Longitude(n)+1.25)*ones(1,26),'k')
  plotm((Latitude(n)-1.25):0.1:(Latitude(n)+1.25),(Longitude(n)-1.25)*ones(1,26),'k')
 end
 setm(gca,'fontweight','bold')
 scatterm(Lat,Lon,50,1:24,'filled','markeredgecolor','k')
 colormap(turbo(numel(Lat)))
 c=colorbar('location','eastoutside');
 caxis([0.5 24.5])
 set(c,'fontweight','bold','ytick',1:24,'yticklabel',[{'San Diego'};{'La Jolla'};{'Los Angeles'};{'Santa Monica'};{'Port San Luis'};{'Monterey'};{'Alameda'};{'San Francisco'};{'Point Reyes'};{'Port Chicago'};{'Arena Cove'};{'Humboldt Bay'};{'Crescent City'};{'Port Orford'};{'Charleston'};{'South Beach'};{'Astoria'};{'Toke Point'};{'Seattle'};{'Port Townsend'};{'Port Angeles'};{'Neah Bay'};{'Friday Harbor'};{'Cherry Point'}])

axes('position',[0.4475 0.725 0.2 0.2])
% origin of your map
originLat = dm2degrees([40 0]);
originLon = dm2degrees([-120 0]);
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','Origin',[originLat originLon]);
geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
plotm((min(Lat)-1):(max(Lat)+1),(min(Lon)-1)*ones(size((min(Lat)-1):(max(Lat)+1))),'k','linewidth',2)
plotm((min(Lat)-1):(max(Lat)+1),(max(Lon)+1)*ones(size((min(Lat)-1):(max(Lat)+1))),'k','linewidth',2)
plotm((min(Lat)-1)*ones(size((min(Lon)-1):(max(Lon)+1))),(min(Lon)-1):(max(Lon)+1),'k','linewidth',2)
plotm((max(Lat)+1)*ones(size((min(Lon)-1):(max(Lon)+1))),(min(Lon)-1):(max(Lon)+1),'k','linewidth',2)
chrisSaveFigurePng(['arsl_f01'],1000)
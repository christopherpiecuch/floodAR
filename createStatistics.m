%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main driver code; performs analyses used to create figures and results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

% define parameters and load NOAA tide gauge IDs
avDa=10; % one-half the 20-day filter cutoff scale
load fileID.mat
ts=datenum(1980,1,1:numel(datenum(1980,1,1):datenum(2016,12,31)));
nt=numel(ts);

% load NOAA tide, water level, and datum
% make daily; remove seasonal and trend
% high-pass filter
for nn=1:numel(ID), disp(['Now on file number ',num2str(nn)])
 % load data, create daily average, keep only 1980-2016
 load(['noaa_tidegauge_',num2str(ID(nn)),'.mat'])
 td(isnan(sl))=nan;
 sl=mean(reshape(sl,24,numel(sl)/24),1); 
 td=mean(reshape(td,24,numel(td)/24),1); 
 dn=floor(mean(reshape(dn,24,numel(dn)/24),1)); 
 sl=sl(find(dn==min(ts)):find(dn==max(ts)));
 td=td(find(dn==min(ts)):find(dn==max(ts)));
 dn=dn(find(dn==min(ts)):find(dn==max(ts)));

 % grab basic stats
 data(nn).na=datum.name';
 data(nn).lo=datum.lon;
 data(nn).la=datum.lat;
 data(nn).tide=td;
 data(nn).date=dn;

 % get sea level on common time grid
 data(nn).sl=sl-td;
 clear dn sl td

 % remove seasonal and trend then high-pass filter (remove 20-day moving median)
 ii=[];
 ii=find(~isnan(data(nn).sl));
 ye=fitharmon_err(ii,data(nn).sl(ii),[365.25 365.25/2],0,1);
 data(nn).sl_dt_ds=data(nn).sl;
 data(nn).sl_dt_ds(ii)=data(nn).sl_dt_ds(ii)-ye';
 data(nn).sl_dt_ds_hf=nan*data(nn).sl_dt_ds;
 data(nn).sl_dt_ds_hf2=nan*data(nn).sl_dt_ds;
 data(nn).sl_dt_ds_hf=data(nn).sl_dt_ds-movmedian(data(nn).sl_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

end

% sort by latitude;
for nn=1:numel(data);
 L(nn)=data(nn).la;
end
[yy ii]=sort(L);
data=data(ii);
clearvars -except TAUX PRECIP avDa data ts nt ID

% load AR catalogue
load('SIO_R1_1948-2017_Comprehensive.mat')
LonLatPairs=unique(Longitude+i*Latitude);

% load slp from ncep-r1 and compute ib colocated with ar catalogue
baseDir='ncepncarr1/';
IB=[];
TAUY=[];
TAUX=[];
PRECIP=[];
% slp
for nn=min(unique(str2num(datestr(ts,10)))):max(unique(str2num(datestr(ts,10)))), disp(num2str(nn))
 % load field
 xx=[]; xx=openNetCDF([baseDir,'slp.',num2str(nn),'.nc']);
 if nn==min(unique(str2num(datestr(ts,10))))
  lo=[]; la=[];
  lo=double(xx.lon);
  la=double(xx.lat);
  % get coordinates
  % initialize arrays
  ib=[];
  for mm=1:numel(LonLatPairs);
   i=[]; i=find(lo==real(LonLatPairs(mm))); I(mm)=i;
   j=[]; j=find(la==imag(LonLatPairs(mm))); J(mm)=j;
   ib(mm,:)=-double(squeeze(xx.slp(I(mm),J(mm),:)))/1029/9.81;
  end
  IB=ib;
  clear lo la ib 
 else
  % initialize arrays
  ib=[];
  for mm=1:numel(LonLatPairs);
   ib(mm,:)=-double(squeeze(xx.slp(I(mm),J(mm),:)))/1029/9.81; % compute local IB response; ignore global-ocean-average slp
  end
  IB=[IB ib];
  clear ib
 end
 clear xx
end

% tauy
for nn=min(unique(str2num(datestr(ts,10)))):max(unique(str2num(datestr(ts,10)))), disp(num2str(nn))
 % load field
 xx=[]; xx=openNetCDF([baseDir,'vflx.sfc.gauss.',num2str(nn),'.nc']);
 if nn==min(unique(str2num(datestr(ts,10))))
  lo=[]; la=[];
  lo=double(xx.lon);
  la=double(xx.lat);
  % get coordinates
  % initialize arrays
  tauv=[];
  for mm=1:numel(LonLatPairs);
   i=[]; i=find(abs(lo-real(LonLatPairs(mm)))==min(abs(lo-real(LonLatPairs(mm))))); I(mm)=i;
   j=[]; j=find(abs(la-imag(LonLatPairs(mm)))==min(abs(la-imag(LonLatPairs(mm))))); J(mm)=j;
   tauy(mm,:)=-double(squeeze(xx.vflx(I(mm),J(mm),:))); % minus sign b/c of weird sign convention in ncep/ncar
  end
  TAUY=tauy;
  clear lo la tauy 
 else
  % initialize arrays
  tauy=[];
  for mm=1:numel(LonLatPairs);
   tauy(mm,:)=-double(squeeze(xx.vflx(I(mm),J(mm),:))); % minus sign b/c of weird sign convention in ncep/ncar
  end
  TAUY=[TAUY tauy];
  clear tauy
 end
 clear xx
end

% taux
for nn=min(unique(str2num(datestr(ts,10)))):max(unique(str2num(datestr(ts,10)))), disp(num2str(nn))
 % load field
 xx=[]; xx=openNetCDF([baseDir,'uflx.sfc.gauss.',num2str(nn),'.nc']);
 if nn==min(unique(str2num(datestr(ts,10))))
  lo=[]; la=[];
  lo=double(xx.lon);
  la=double(xx.lat);
  % get coordinates
  % initialize arrays
  tauv=[];
  for mm=1:numel(LonLatPairs);
   i=[]; i=find(abs(lo-real(LonLatPairs(mm)))==min(abs(lo-real(LonLatPairs(mm))))); I(mm)=i;
   j=[]; j=find(abs(la-imag(LonLatPairs(mm)))==min(abs(la-imag(LonLatPairs(mm))))); J(mm)=j;
   taux(mm,:)=-double(squeeze(xx.uflx(I(mm),J(mm),:))); % minus sign b/c of weird sign convention in ncep/ncar
  end
  TAUX=taux;
  clear lo la taux 
 else
  % initialize arrays
  taux=[];
  for mm=1:numel(LonLatPairs);
   taux(mm,:)=-double(squeeze(xx.uflx(I(mm),J(mm),:))); % minus sign b/c of weird sign convention in ncep/ncar
  end
  TAUX=[TAUX taux];
  clear taux
 end
 clear xx
end

% precip
for nn=min(unique(str2num(datestr(ts,10)))):max(unique(str2num(datestr(ts,10)))), disp(num2str(nn))
 % load field
 xx=[]; xx=openNetCDF([baseDir,'prate.sfc.gauss.',num2str(nn),'.nc']);
 if nn==min(unique(str2num(datestr(ts,10))))
  lo=[]; la=[];
  lo=double(xx.lon);
  la=double(xx.lat);
  % get coordinates
  % initialize arrays
  prec=[];
  for mm=1:numel(LonLatPairs);
   i=[]; i=find(abs(lo-real(LonLatPairs(mm)))==min(abs(lo-real(LonLatPairs(mm))))); I(mm)=i;
   j=[]; j=find(abs(la-imag(LonLatPairs(mm)))==min(abs(la-imag(LonLatPairs(mm))))); J(mm)=j;
   precip(mm,:)=double(squeeze(xx.prate(I(mm),J(mm),:))); 
  end
  PRECIP=precip;
  clear lo la precip 
 else
  % initialize arrays
  precip=[];
  for mm=1:numel(LonLatPairs);
   precip(mm,:)=double(squeeze(xx.prate(I(mm),J(mm),:))); 
  end
  PRECIP=[PRECIP precip];
  clear precip
 end
 clear xx

end

% remove time-mean values from everything
for n=1:numel(LonLatPairs), disp(num2str(n))
 IB(n,:)=detrend(IB(n,:),'constant');
 TAUY(n,:)=detrend(TAUY(n,:),'constant');
 TAUX(n,:)=detrend(TAUX(n,:),'constant');
 PRECIP(n,:)=detrend(PRECIP(n,:),'constant');
end
clearvars -except TAUX PRECIP avDa data IB TAUY LonLatPairs ts nt ID

% reload AR catalogue
load('SIO_R1_1948-2017_Comprehensive.mat')
% as in shinoda et al., exclude all events for which IVT<500 kg m−1 s−1
%ii=find(IVT<500);
ii=find(IVT<500|Year<1980|Year>2016);
Day(ii)=[];
Hour(ii)=[];
IVT(ii)=[];
IWT(ii)=[];
Latitude(ii)=[];
Longitude(ii)=[];
Month(ii)=[];
Number(ii)=[];
Uwind(ii)=[];
Vwind(ii)=[];
Year(ii)=[];
clear ii
LonLatPairs=LonLatPairs-360;
Longitude=Longitude-360;

% loop back through the tide gauge data and find closest ncep-r1 grid cell
for nn=1:numel(data), disp(num2str(nn))
tic
 d=[]; d=distance(data(nn).la,data(nn).lo,imag(LonLatPairs),real(LonLatPairs),6371);
 i=[]; i=find(d==min(d),1,'first');
 j=[]; j=find(Longitude==real(LonLatPairs(i))&Latitude==imag(LonLatPairs(i)));

 % grab AR & IB and filter just like real data
 data(nn).i=i;
 data(nn).ar_da=Day(j);
 data(nn).ar_hr=Hour(j);
 data(nn).ar_mo=Month(j);
 data(nn).ar_yr=Year(j);
 data(nn).ar_no=Number(j);
 data(nn).ar_uw=Uwind(j);
 data(nn).ar_vw=Vwind(j);
 data(nn).ar_iv=IVT(j);
 data(nn).ar_iw=IWT(j);
 data(nn).ib=squeeze(IB(i,:))';
 data(nn).ib(find(isnan(data(nn).sl)))=nan;
 data(nn).tauy=squeeze(TAUY(i,:))';
 data(nn).tauy(find(isnan(data(nn).sl)))=nan;
 data(nn).taux=squeeze(TAUX(i,:))';
 data(nn).taux(find(isnan(data(nn).sl)))=nan;
 data(nn).precip=squeeze(PRECIP(i,:))';
 data(nn).precip(find(isnan(data(nn).sl)))=nan;

 % compute hilbert transforms
 data(nn).Hib=imag(hilbert(squeeze(IB(i,:))'));
 data(nn).Hib(find(isnan(data(nn).sl)))=nan;
 data(nn).Htauy=imag(hilbert(squeeze(TAUY(i,:))'));
 data(nn).Htauy(find(isnan(data(nn).sl)))=nan;
 data(nn).Htaux=imag(hilbert(squeeze(TAUX(i,:))'));
 data(nn).Htaux(find(isnan(data(nn).sl)))=nan;
 data(nn).Hprecip=imag(hilbert(squeeze(PRECIP(i,:))'));
 data(nn).Hprecip(find(isnan(data(nn).sl)))=nan;

 % do filtering; remove seasonal and trend then high pass
 ii=[];
 ii=find(~isnan(data(nn).ib));
 ye=fitharmon_err(ii,data(nn).ib(ii),[365.25 365.25/2],0,1);
 data(nn).ib_dt_ds=data(nn).ib;
 data(nn).ib_dt_ds(ii)=data(nn).ib_dt_ds(ii)-ye;
 data(nn).ib_dt_ds_hf=nan*data(nn).ib_dt_ds;
 data(nn).ib_dt_ds_hf=data(nn).ib_dt_ds-movmedian(data(nn).ib_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).tauy));
 ye=fitharmon_err(ii,data(nn).tauy(ii),[365.25 365.25/2],0,1);
 data(nn).tauy_dt_ds=data(nn).tauy;
 data(nn).tauy_dt_ds(ii)=data(nn).tauy_dt_ds(ii)-ye;
 data(nn).tauy_dt_ds_hf=nan*data(nn).tauy_dt_ds;
 data(nn).tauy_dt_ds_hf=data(nn).tauy_dt_ds-movmedian(data(nn).tauy_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).taux));
 ye=fitharmon_err(ii,data(nn).taux(ii),[365.25 365.25/2],0,1);
 data(nn).taux_dt_ds=data(nn).taux;
 data(nn).taux_dt_ds(ii)=data(nn).taux_dt_ds(ii)-ye;
 data(nn).taux_dt_ds_hf=nan*data(nn).taux_dt_ds;
 data(nn).taux_dt_ds_hf=data(nn).taux_dt_ds-movmedian(data(nn).taux_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).precip));
 ye=fitharmon_err(ii,data(nn).precip(ii),[365.25 365.25/2],0,1);
 data(nn).precip_dt_ds=data(nn).precip;
 data(nn).precip_dt_ds(ii)=data(nn).precip_dt_ds(ii)-ye;
 data(nn).precip_dt_ds_hf=nan*data(nn).precip_dt_ds;
 data(nn).precip_dt_ds_hf=data(nn).precip_dt_ds-movmedian(data(nn).precip_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).Hib));
 ye=fitharmon_err(ii,data(nn).Hib(ii),[365.25 365.25/2],0,1);
 data(nn).Hib_dt_ds=data(nn).Hib;
 data(nn).Hib_dt_ds(ii)=data(nn).Hib_dt_ds(ii)-ye;
 data(nn).Hib_dt_ds_hf=nan*data(nn).Hib_dt_ds;
 data(nn).Hib_dt_ds_hf=data(nn).Hib_dt_ds-movmedian(data(nn).Hib_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).Htauy));
 ye=fitharmon_err(ii,data(nn).Htauy(ii),[365.25 365.25/2],0,1);
 data(nn).Htauy_dt_ds=data(nn).Htauy;
 data(nn).Htauy_dt_ds(ii)=data(nn).Htauy_dt_ds(ii)-ye;
 data(nn).Htauy_dt_ds_hf=nan*data(nn).Htauy_dt_ds;
 data(nn).Htauy_dt_ds_hf=data(nn).Htauy_dt_ds-movmedian(data(nn).Htauy_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).Htaux));
 ye=fitharmon_err(ii,data(nn).Htaux(ii),[365.25 365.25/2],0,1);
 data(nn).Htaux_dt_ds=data(nn).Htaux;
 data(nn).Htaux_dt_ds(ii)=data(nn).Htaux_dt_ds(ii)-ye;
 data(nn).Htaux_dt_ds_hf=nan*data(nn).Htaux_dt_ds;
 data(nn).Htaux_dt_ds_hf=data(nn).Htaux_dt_ds-movmedian(data(nn).Htaux_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

 ii=[];
 ii=find(~isnan(data(nn).Hprecip));
 ye=fitharmon_err(ii,data(nn).Hprecip(ii),[365.25 365.25/2],0,1);
 data(nn).Hprecip_dt_ds=data(nn).Hprecip;
 data(nn).Hprecip_dt_ds(ii)=data(nn).Hprecip_dt_ds(ii)-ye;
 data(nn).Hprecip_dt_ds_hf=nan*data(nn).Hprecip_dt_ds;
 data(nn).Hprecip_dt_ds_hf=data(nn).Hprecip_dt_ds-movmedian(data(nn).Hprecip_dt_ds,2*avDa,'omitnan','Endpoints','shrink');

toc
end
%clearvars -except TAUX PRECIP avDa data

% now, for each TG site and AR number, find the day when the AR was peak by IVT  
for nn=1:numel(data), disp(num2str(nn))
tic
 uq=[]; uq=unique(data(nn).ar_no);
 for mm=1:numel(uq)
  data(nn).lengthscale(mm)=distance(data(nn).la,data(nn).lo,min(Latitude(find(Number==uq(mm)))),data(nn).lo,6371);
  cond1=[]; cond2=[];
  cond1=data(nn).ar_no==uq(mm);
  cond2=data(nn).ar_iv==max(data(nn).ar_iv(find(data(nn).ar_no==uq(mm))));
  ii=[]; ii=find(cond1&cond2,1,'first');
  data(nn).ar_uq_iw(mm)=data(nn).ar_iw(ii);
  data(nn).ar_uq_iv(mm)=data(nn).ar_iv(ii);
  data(nn).ar_uq_uw(mm)=data(nn).ar_uw(ii);
  data(nn).ar_uq_vw(mm)=data(nn).ar_vw(ii);

  la1=[]; lo1=[];
  la1=data(nn).la;
  la2=min(Latitude(find(Number==uq(mm))));
  lo2=lo1; % same as first; only wanting latitudinal distance

  ll=[]; ll=Latitude(find(Number==uq(mm)));
  vv=[]; vv=IVT(find(Number==uq(mm)));
  la3=ll(find(vv==max(vv)));

  data(nn).ar_dist2first(mm)=la1-la2;
  data(nn).ar_dist2max(mm)=la1-la3;
  data(nn).ar_uq_no(mm)=uq(mm);
  data(nn).ar_uq_dn(mm)=datenum(data(nn).ar_yr(ii),data(nn).ar_mo(ii),data(nn).ar_da(ii));
  jj=find(data(nn).ar_uq_dn(mm)==ts);
  data(nn).ar_uq_tg_ind(mm)=jj;
  % consider points within +/-5 days of AR landfall
  if jj>(1+9)&jj<(numel(datenum(1980,1,1):datenum(2016,12,31))-9)
   data(nn).ar_uq_sl_window(mm,:)=data(nn).sl_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_ib_window(mm,:)=data(nn).ib_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_tauy_window(mm,:)=data(nn).tauy_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_taux_window(mm,:)=data(nn).taux_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_precip_window(mm,:)=data(nn).precip_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_Hib_window(mm,:)=data(nn).Hib_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_Htauy_window(mm,:)=data(nn).Htauy_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_Htaux_window(mm,:)=data(nn).Htaux_dt_ds_hf((jj-5):(jj+5));
   data(nn).ar_uq_Hprecip_window(mm,:)=data(nn).Hprecip_dt_ds_hf((jj-5):(jj+5));
  else
   data(nn).ar_uq_sl_window(mm,:)=nan(11,1);
   data(nn).ar_uq_ib_window(mm,:)=nan(11,1);
   data(nn).ar_uq_tauy_window(mm,:)=nan(11,1);
   data(nn).ar_uq_taux_window(mm,:)=nan(11,1);
   data(nn).ar_uq_precip_window(mm,:)=nan(11,1);
   data(nn).ar_uq_sl_window(mm,:)=nan(11,1);
   data(nn).ar_uq_Hib_window(mm,:)=nan(11,1);
   data(nn).ar_uq_Htauy_window(mm,:)=nan(11,1);
   data(nn).ar_uq_Htaux_window(mm,:)=nan(11,1);
   data(nn).ar_uq_Hprecip_window(mm,:)=nan(11,1);
  end
 end
 kk=[]; kk=find(isnan(data(nn).ar_uq_sl_window(:,6)));
 data(nn).ar_uq_sl_window(kk,:)=[];
 data(nn).ar_uq_ib_window(kk,:)=[];
 data(nn).ar_uq_tauy_window(kk,:)=[];
 data(nn).ar_uq_taux_window(kk,:)=[];
 data(nn).ar_uq_precip_window(kk,:)=[];
 data(nn).ar_uq_Hib_window(kk,:)=[];
 data(nn).ar_uq_Htauy_window(kk,:)=[];
 data(nn).ar_uq_Htaux_window(kk,:)=[];
 data(nn).ar_uq_Hprecip_window(kk,:)=[];
 data(nn).lengthscale(kk)=[];
 data(nn).ar_uq_tg_ind(kk)=[];
 data(nn).ar_uq_dn(kk)=[];
 data(nn).ar_uq_no(kk)=[];

 % ok, what did you do? you picked the number of unique atmospheric rivers such that at their peak you have sea level data
 data(nn).ar_dist2first(kk)=[];
 data(nn).ar_uq_iw(kk)=[];
 data(nn).ar_uq_iv(kk)=[];
 data(nn).ar_uq_uw(kk)=[];
 data(nn).ar_uq_vw(kk)=[];
toc
end
clearvars -except TAUX PRECIP avDa data ts

for nn=1:numel(data)
 Lo(nn)=data(nn).lo;
 La(nn)=data(nn).la;
 Md(nn)=mean(data(nn).ar_dist2first);
 Mm(nn)=mean(data(nn).ar_dist2max);
end

rng(123);
% take values from day contemporaneous with AR landfall
% (6th element of 11-element +/-5 day array)
for nn=1:numel(data), disp(num2str(nn))
 results(nn).data=data(nn).ar_uq_sl_window(:,6);
 results(nn).tauy=data(nn).ar_uq_tauy_window(:,6);
 results(nn).taux=data(nn).ar_uq_taux_window(:,6);
 results(nn).pres=data(nn).ar_uq_ib_window(:,6);
 results(nn).precip=data(nn).ar_uq_precip_window(:,6);
 results(nn).Htauy=data(nn).ar_uq_Htauy_window(:,6);
 results(nn).Htaux=data(nn).ar_uq_Htaux_window(:,6);
 results(nn).Hpres=data(nn).ar_uq_Hib_window(:,6);
 results(nn).Hprecip=data(nn).ar_uq_Hprecip_window(:,6);
 results(nn).data_max=max(data(nn).ar_uq_sl_window');
 results(nn).tauy_max=max(data(nn).ar_uq_tauy_window');
 results(nn).taux_max=max(data(nn).ar_uq_taux_window');
 results(nn).pres_max=max(data(nn).ar_uq_ib_window');
 results(nn).precip_max=max(data(nn).ar_uq_precip_window');
end

for nn=1:numel(data)
 arss(nn).name=data(nn).na';
 arss(nn).longitude=data(nn).lo;
 arss(nn).latitude=data(nn).la;
 arss(nn).data=results(nn).data;
 arss(nn).precip=results(nn).precip/1029; % make precip in m/s units
 arss(nn).tauy=results(nn).tauy;
 arss(nn).taux=results(nn).taux;
 arss(nn).pres=results(nn).pres*(-1029*9.81); % convert back to hPa
 arss(nn).Hprecip=results(nn).Hprecip/1029;
 arss(nn).Htauy=results(nn).Htauy;
 arss(nn).Htaux=results(nn).Htaux;
 arss(nn).Hpres=results(nn).Hpres*(-1029*9.81);
 arss(nn).data_max=results(nn).data_max;
 arss(nn).precip_max=results(nn).precip_max/1029;
 arss(nn).tauy_max=results(nn).tauy_max;
 arss(nn).taux_max=results(nn).taux_max;
 arss(nn).pres_max=results(nn).pres_max*(-1029*9.81);
end

% save out
save(['ar_and_ss_statistics_w_precip_',num2str(avDa),'.mat'],'arss')
save(['sealevel_statistics_w_precip_',num2str(avDa),'.mat'],'data')
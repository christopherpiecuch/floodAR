%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for downloading tide-gauge station datums, sea-level observations,
% and tidal preditions from NOAA NOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
% NOAA IDs of tide gauges used in this study (cf. Table S1)
ID=[9410170 9410230 9410660 9410840 9412110 9413450 9414290 9414750 9415020 9415144 9416841 9418767 9431647 9432780 9435380 9440910 9444090 9444900 9447130 9449424 9449880 9419750 9439040 9443090];
save(['fileID.mat'],'ID')
t0=1980;
tf=2016;

% loop through, load data, save
for id=1:numel(ID), %disp(num2str(id))
 tic
 [datum]=noaaDatums(ID(id));
 [sl,dn]=noaaSealevel(ID(id),t0,tf); % hourly observed water levels
 [td,dn]=noaaTide(ID(id),t0,tf); % hourly observed water levels
 save(['noaa_tidegauge_',num2str(datum.id),'.mat'],'datum','sl','td','dn')
 clear datum sl td dn
 toc
end
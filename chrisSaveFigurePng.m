%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quick code to save out high-resolution PNG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chrisSaveFigurePng(fileName,rez)
%rez=resolution (dpi) of final graphic
f=gcf; %f is the handle of the figure you want to export
figpos=getpixelposition(f); 
resolution=get(0,'ScreenPixelsPerInch'); 
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fileName,'-dpng',['-r',num2str(rez)],'-painters') %save file 
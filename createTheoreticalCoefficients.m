%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of theoretical regression coefficients (as and bs) shown in
% Figure S2 and described in Supporting Information Text S4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

% define constants; cf. Table S2
rho=1000; % kg/m^3 density
g=10; % m/s^2 gravity
H=[100 200]; % m depth of shelf
f=[0.6e-4 1.1e-4]; % 1/s typical midlatitude Coriolis parameter value
% define uncertain parameters
k=[1/(50*1e3) 1/(200*1e3)]; % 1/m short and long decay scales
R=[1e-4 1e-2]; % m/s strong and weak friction

% set timescales
Tmin=1*86400; % 1 days minimum
Tmax=6*86400; % 6 days maximum
sigmaMax=2*pi/Tmin;
sigmaMin=2*pi/Tmax;
sigma=sigmaMin:(sigmaMin/10):sigmaMax;

for ns=1:numel(sigma), %disp(num2str(ns))
 for nk=1:numel(k)
  for nr=1:numel(R)
   for nh=1:numel(H)
    for nf=1:numel(f)
	lambda=[]; lambda=R(nr)/H(nh); % 1/s inverse frictional timescale
    LR=[]; LR=sqrt(g*H(nh))/f(nf); % m Rossby radius
	r=[]; r=(1+((lambda)/(sigma(ns)))^2)^(-1/4);
	theta=[]; theta=0.5*mod(atan2(-lambda/sigma(ns),1),2*pi);
	z=[]; z=r*exp(i*theta);
	kappa=[]; kappa=z/LR;
	denom=[]; denom=(1/((k(nk)-kappa)));

    % compute coefficients; see Equations S17-S20
	% x wind
	xx=[]; xx=denom*(1/(rho*g*H(nh)));
	ax(ns,nk,nr,nh,nf)=real(xx);
	bx(ns,nk,nr,nh,nf)=imag(xx);

	% y wind
	xx=[]; xx=denom*(f(nf)/(rho*g*H(nh)))*((lambda+i*sigma(ns))/(lambda^2+sigma(ns)^2));
	ay(ns,nk,nr,nh,nf)=real(xx);
	by(ns,nk,nr,nh,nf)=imag(xx);

	% pressure
	xx=[]; xx=denom*(-(k(nk))/(rho*g));
	ap(ns,nk,nr,nh,nf)=real(xx);
	bp(ns,nk,nr,nh,nf)=imag(xx);

	% rainfall
	xx=[]; xx=-i*denom*(kappa/sigma(ns));
	aq(ns,nk,nr,nh,nf)=real(xx);
	bq(ns,nk,nr,nh,nf)=imag(xx);

    end
   end
  end % for nr=1:numel(lambda)
 end % for nk=1:numel(k)
end % for ns=1:numel(sigma)

% reshape and take maximum and minumum values
ax=reshape(ax,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
bx=reshape(bx,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
ay=reshape(ay,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
by=reshape(by,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
ap=reshape(ap,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
bp=reshape(bp,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
aq=reshape(aq,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
bq=reshape(bq,numel(ap)/(numel(f)*numel(k)*numel(H)*numel(R)),(numel(f)*numel(k)*numel(H)*numel(R)));
ap=[min(ap(:)) max(ap(:))];
bp=[min(bp(:)) max(bp(:))];
aq=[min(aq(:)) max(aq(:))];
bq=[min(bq(:)) max(bq(:))];
ax=[min(ax(:)) max(ax(:))];
bx=[min(bx(:)) max(bx(:))];
ay=[min(ay(:)) max(ay(:))];
by=[min(by(:)) max(by(:))];

clear ans
save('theoretical_coefficients','a*','b*')
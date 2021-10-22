%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piecuch, C. G., et al. (2021)
% High-Tide Floods and Storm Surges During Atmospheric Rivers on the US West Coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simultaneous OLS fit of linear trend plus sinuoids.
% Used for removing long-term rate of change and seasonal cycle.
% Original code written by Dr. Katherine J. Quinn 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ye,amp,amp_err,phase,phase_err,yint,yint_err,slope,slope_err] = fitharmon_err(x,y,period,x0,linfit,ysig)
% [ye,amp,amp_err,phase,phase_err,yint,yint_err,slope,slope_err] = fitharmon_err(x,y,period,x0,linfit,ysig)
% fit harmonic function
% period can be vector of periods to fit to
% i.e. annual, semiannual, quarter-annual.
% NB: period units will be same as x units.
% x0 is the point to calculate phase from
% linfit is flag to fit for trend
% 0 = none, 1 = trend
% Solution in form:  Asin(2*pi/period*(x-x0) + phase*pi/180)
% errors are standard errors
% 95% CI = stderr*1.96

% Output:
% ye = y estimated (fitted solution)
% amp = amplitudes
% phase = phases (in degrees)
% slope = slope of linear fit
np = length(period);
x = x(:);
y = y(:);
my = mean(y);
y = y - my;

if nargin < 6
  ysig = ones(size(y));
end
ysig = ysig(:);

% set up jacobian
A = [];
for i = 1:np
  A = [A, sin(2*pi./period(i).*(x - x0)), cos(2*pi./period(i).*(x - x0))];
end
if linfit == 1
  A = [A, ones(size(x)), x];
else
  A = [A, ones(size(x))];
end

% solve linear least squares: y = A*b
[b,b_err] = lscov(A,y,ysig);

ye = A*b + my;
yint = b(2*np+1) + my;
yint_err = b_err(2*np+1);
if linfit == 1
  slope = b(2*np+2);
  slope_err = b_err(2*np+2);
end
amp = []; amp_err = [];
phase = []; phase_err = [];
for i = 1:np
  c = b(2*i-1);
  s = b(2*i);
  a = sqrt(c*c + s*s);
  p = atan2(s/a,c/a)*180/pi;
  if p < 0, p = p+360; end
  amp = [amp; a];
  phase = [phase; p];
  c_err = b_err(2*i-1);
  s_err = b_err(2*i);
  a_err = sqrt( (s/a*s_err).^2 + (c/a*c_err).^2 );
  p_err = sqrt( (c/a/a*s_err).^2 + (s/a/a*c_err).^2 ) *180/pi;
  amp_err = [amp_err; a_err];
  phase_err = [phase_err; p_err];
end
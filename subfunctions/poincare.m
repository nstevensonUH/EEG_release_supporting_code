function [sd1, sd2] = poincare(rr);
% function [sd1, sd2] = poincare(rr, msd);
%
% INPUT: rr - is the interpolated NN interval
%        msd - missing data vector 
%
% OUTPUT: sd1 -  
%         sd2 -
%
% Nathan Stevenson
% Neonatal Brain Research Group
% 15/02/2013

x = rr(1:end-3)*1000;
y = rr(4:end)*1000;
xdash = x.*cos(pi/4)+y.*sin(pi/4);
ydash = x.*sin(pi/4)-y.*cos(pi/4);
sd1 = 4*std(ydash); sd2 = 4*std(xdash);

function [o,v]=tdoa_offset_one(utmp0,u_anchor0,Abase0);
%[o,v]=tdoa_offset_one(utmp0,u_anchor0,Abase0);
%

% blubb0 = (utmp0-oo).^2 - (u_anchor0-o_anchor).^2;
% but since u_anchor is already u_anchor0-o_anchor this becomes
% blubb0 = (utmp0-oo).^2 - u_anchor0.^2;
% = utmp0.^2 - 2*utmp0*oo - u_anchor0.^2;


blubb0 = utmp0.^2 - u_anchor0.^2;
blubb0 = blubb0(2:end)-blubb0(1);
blubb1 = -2*utmp0;
blubb1 = blubb1(2:end)-blubb1(1);
% Ok Now we should solve 
% blubb0+blubb1*oo = Abase0*vv
% or 
% [Abase blubb1]*[vv;-oo]=blubb0;

tmp = [Abase0 blubb1]\blubb0;
v = tmp(1:size(Abase0,2));
o = -tmp(end);



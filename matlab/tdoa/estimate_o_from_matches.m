function [matches,res] = estiamate_o_from_matches(rawmatches,settings);
%[matches,res] = estiamate_o_from_matches(rawmatches,settings);
%

u = rawmatches.u;
ok = rawmatches.uok;
[o,v,uuu,A,inliers_ut,pp,cc]=tdoa_offset_ransac(u,ok);

% % Check nr of inliers????
ok = ( (sum(inliers_ut)>5) & (o<0) );
% Also remove outliers in terms of o and A
blubb = sort(o(find(ok)));
Tlo = blubb( round(length(blubb)*0.1));
Thi = blubb( round(length(blubb)*0.9));
ok = ( (sum(inliers_ut)>5) & (o<Thi) & (o>Tlo) & (o<min(u)) );

%jsel = find(sum(inliers_ut)>5);  % Also remove outliers in terms of o and A
jsel = find( ok );  % Also remove outliers in terms of o and A

Atmp = A(:,jsel);
[minv,mini]=min(sum(abs(Atmp)));
janchor = jsel(mini);
jsel(mini)=[];
jsel = [janchor jsel];

% Extract the part of the data (uij) that has been matched using ransac
utmp = rawmatches.u(:,jsel);
oktmp = rawmatches.uok(:,jsel);
Atmp = A(:,jsel);
Atmp = Atmp(:,2:end);
utimes = rawmatches.utimes(jsel);
uindex = rawmatches.uindex(jsel);
%vtmp = v(:,jsel);
ojoj = find(isnan(utmp));
utmp(ojoj)=zeros(size(ojoj));
itmp = inliers_ut(:,jsel);
otmp = o(jsel);

if 0, % Checking for errors
    [d,dc,d2,d2c,d2r,d2rc]=checkequations(utmp,itmp,otmp,Atmp);    
end

%%hack%%%
oopt = otmp;
o = otmp;
Aopt = Atmp;
[paramopt.m,paramopt.n] = size(utmp);
paramopt.D = utmp;
resopt = '?';
jacopt = '?';
%%%%%%%%%
% [oopt,Aopt,paramopt,resopt,jacopt]=bundle_rank1(utmp,itmp,otmp',3,Atmp);
% norm(resopt)
% %[oopt,Aopt,paramopt,resopt,jacopt]=bundle_rank1(utmp,itmp,oopt+1,3,Aopt);
% %norm(resopt)
% [oopt,Aopt,paramopt,resopt,jacopt]=bundle_rank1(utmp,itmp,oopt,3,Aopt);
% norm(resopt)

% o = oopt';

if 0, % Checking for errors
    [d,dc,d2,d2c,d2r,d2rc]=checkequations(utmp,itmp,oopt',Aopt);    
end


matches.u = utmp;
matches.uok = oktmp;
matches.uindex = uindex;
matches.utimes = utimes;
matches.uinliers = itmp;
res.A = Aopt;
res.jsel = jsel;
res.janchor = janchor;
res.o = oopt;
res.paramopt = paramopt;
res.resopt = resopt;
res.jacopt = jacopt;
res.pp = pp;
res.cc = cc;



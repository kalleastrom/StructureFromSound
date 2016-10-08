function [data]=sampleWithU(a,t0,uvc,w,settings);
%
%


t0i = t0*settings.sr;
uvci = uvc*settings.sr/settings.v;

t0i = round(t0i);
uvci = round(uvci);

data = zeros(size(a,1),2*w+1);

for kk = 1:size(a,1);
    data(kk,:)=a(kk,(t0i+uvci(kk)-w):(t0i+uvci(kk)+w));
end;


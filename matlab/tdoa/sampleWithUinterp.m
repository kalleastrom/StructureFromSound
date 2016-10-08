function [data,datap]=sampleWithUinterp(a,t0,uvc,w,settings);
%
%

t0i = t0*settings.sr; %sample index
uvci = uvc*settings.sr/settings.v; %matches(by pixel)

%t0i = round(t0i);
%uvci = round(uvci);

data0 = zeros(size(a,1),2*w+1);
data = zeros(size(a,1),2*w+1);
for kk = 1:size(a,1);
    theshift = t0i+uvci(kk);
    theshift_integerpart = round(theshift);
    theshift_rest = theshift - theshift_integerpart;
    
    %data0(kk,:)=a(kk,(theshift_integerpart-w):(theshift_integerpart+w)); 
    
    tmp1 = a(kk,(theshift_integerpart-w):(theshift_integerpart+w));
  

    x0 = [-5:5];
    sigma = 1.2;
    miu = -theshift_rest; % or maybe -theshift_rest
    h = (1/sqrt(2*pi)*sigma)*exp(-(x0 - miu).^2/(2*sigma^2));
    h = h/sum(h);
    hprime = (-2*(x0-miu)/(2*sigma^2)).*h;
    
    tmp2 = conv(tmp1,h,'same');
    tmp2prime = conv(tmp1,hprime,'same');
    
    data(kk,:) = tmp2;
    datap(kk,:) = tmp2prime;
end;

datap = datap/(settings.v/settings.sr);

%datap = [diff(data')' zeros(8,1)]/(settings.v/settings.sr);
% The row above is almost correct. 



    % The value t0i+uvci(kk) is not always integer. 
    % perhaps take the nearest integer value first here, but 
    %
    
    %Ideal interpolation
    %x = [1:size(data0,2)];
    %y = data0(kk,:);
    %xi = [1:.5:size(data0,2)];
    %yi = interp1(x,y,xi);
    %figure(3);clf
    %plot(x,y,'o',xi,yi)
    
    %Smoothing with gaussian
%     x0 = [-5:5];
%     sigma = 1.2;
%     miu = 0.2;
%     h = (1/sqrt(2*pi)*sigma)*exp(-(x0 - miu).^2/(2*sigma^2));
%     h = h/sum(h);
%     h0 = (1/sqrt(2*pi)*sigma)*exp(-(x0 - 0).^2/(2*sigma^2));
%     h0 = h0/sum(h0);
%     figure(5),
%     hold off;
%     plot(x0,h);
%     hold on;
%     plot(x0,h0,'r');
%     
%     tmp2 = conv(tmp1,h,'same');
%     tmp0 = conv(tmp1,h0,'same');
%     figure(6);
%     hold off;
%     plot(tmp1,'b');
%     hold on;
%     plot(tmp2,'g');
%     plot(tmp0,'r');
%     [tmp1;tmp2;tmp0]

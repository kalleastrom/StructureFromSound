function [cres] = myplot(rns,cmatches,cres,settings);

x = cres.x;
y = cres.y;
o = cres.o;

utimes = cmatches.utimes;
u = cmatches.u;

dd = reshape(sqrt( sum( (kron(ones(1,size(y,2)),x) - ...
    kron(y,ones(1,size(x,2))) ).^2 , 1 ) ),size(x,2),size(y,2)) + repmat(o,size(x,2),1);

[I,J]=find(cmatches.uinliers);
D = u(find(cmatches.uinliers));
[res,Dproj]=calcres(D,I,J,x,y,o);
resm = full(sparse(I,J,res));
cres.res = res;
cres.resm = resm;

[maxv,maxi]=max(sum(abs(resm)));

if settings.doplot,
    figure(3);
    plot(res*1000,'b.');
    figure(4); clf;
    hist(res(:)*settings.sr/settings.v,100);
    title(['Histogram of errors in terms of sampling points']);
    figure(5); clf;
    hist(res(:)*1000,100);
    title(['Histogram of errors in millimeter']);
end;

if settings.doplot,
    for kk = 2:size(cmatches.u,1);
        figure(10+kk);
        clf;
        colormap(gray);
        imagesc(settings.xtid,settings.ymeter,rns{kk});
        xlabel('time in seconds');
        ylabel('offset in meters');
        hold on;
        plot(utimes,dd(kk,:)-dd(1,:),'b*'); % Haven't kept track of xxi
        plot(utimes(maxi),dd(kk,maxi)-dd(1,maxi),'r*');
        plot(utimes,u(kk,:),'go');
        plot(utimes(maxi),u(kk,maxi),'yo');
        %keyboard;
        %pause;
        %drawnow;
    end
end;

if settings.doplot,
    figure(2);
    clf; hold off;
    rita3(cres.y,'r-');
    hold on
    rita3(cres.x,'b*');
    for mi = 1:settings.mm;
        text(cres.x(1,mi),cres.x(2,mi),cres.x(3,mi),num2str(mi),'FontSize',20);
    end
end


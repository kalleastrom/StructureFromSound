%% Define settings

settings.v = 340; % Speed of sound
settings.in = 1000; % nr of search positions (x-axis in rns-plot)
settings.sw = 500; % Search width in sampling points (y-axis in rns-plot)
settings.swstep = 3; % Search steps
settings.wlist = [2000]; % Window size used in matching
settings.scorefunction = 'score1';
settings.xinorm = [1 2 4];
%
settings.doplot = 1;
settings.doverbose = 1;
settings.doprint = 0;

%% Read files

clear a
mm = 8;
%savefolder = '/Users/kalle/Documents/projekt/antennaResection/data/experiment20140918/export/'; expname = 'bassh3-';
%savefolder = '/Users/kalle/Documents/projekt/antennaResection/data/experimentEhuset_20140127/'; expname = 'calib1-';
%savefolder = '/Users/kalle/Documents/projekt/antennaResection/data/experimentEhuset20150205/export/'; expname = 'Musik_mjuka_rorelser-';
%savefolder = '/Users/kalle/Documents/projekt/antennaResection/data/experimentEhuset20150205/export/'; expname = fileNameBase;
%savefolder = 'R:\experiment20140918\export\'; expname = 'bassh3-';
%savefolder = 'R:\experiment20141127\export\'; expname = 'grieg5-'
[a,fs] = readaudio([savefolder expname],'.aiff',mm,1:mm);

%% Calculate additional variables in settings
settings.sr = fs;
settings.nn = size(a,2);
settings.tt = settings.nn/settings.sr;
settings.isel = floor((1:(settings.in-1))*settings.nn/settings.in);
settings.dsel = (-settings.sw):settings.swstep:settings.sw;
settings.xtid = settings.isel/settings.sr;
settings.ymeter = settings.dsel*340/settings.sr;

%% Calculate correlation patterns (PART A)
partasavename = [savefolder expname 'A.mat'];
if ~exist(partasavename),
    [rns,rawmatches] = calculate_matching_score_and_matches_1(a,settings);
    eval(['save ' partasavename ' rns settings rawmatches']);
else
    eval(['load ' partasavename]);
end

% % %hack för bassh3:
% % rawmatches.uok(7,1420:1550) = 0;
% % rawmatches.u(7,1420:1550) = NaN;
% % rawmatches.uij(7,1420:1550) = NaN;


%% plots

if settings.doplot
    tmp = rns{4};
    tmp = tmp./repmat(max(tmp),size(tmp,1),1);
    figure(1);
    clf;
    hold off;
    colormap(gray);
    imagesc(settings.xtid,settings.ymeter,-tmp);
    xlabel('Position, t, along track 1 in seconds');
    ylabel('Offsets, u, in meter');
    title('Matching score for channel 1 vs channel 4');
    printFig('matchingScoreCureChoir_1_4');
end

%% Estimate offsets from rawmatches (PART B)
%settings.in = 1999;
partbsavename = [savefolder expname 'B.mat'];
if ~exist(partbsavename),
    [bmatches,o,bres] = estimate_o_from_matches(rawmatches,settings);
    eval(['save ' partbsavename ' bmatches o bres']);
else
    eval(['load ' partbsavename]);
end

%myplot(rns,bmatches,bres,settings);

%% plots

if 0, %settings.doplot
    tmp1 = rawmatches.uok;
    tmp2 = bmatches.uinliers;
    ui1  = rawmatches.uindex;
    ui2  = bmatches.uindex;
    perm1 = zeros(1,size(ui2,2));
    for kk = 1:size(ui2,2);
        perm1(kk) = find(ui1==ui2(kk));
    end;
    tmp22 = zeros(size(tmp1));
    tmp22(:,perm1)=tmp2;
    tmp3 = zeros(size(tmp1));
    tmp3(bres.cc,bres.pp)=ones(size(7,6));
    [i1,j1]=find(tmp1-tmp22);
    [i2,j2]=find(tmp22-tmp3);
    [i3,j3]=find(tmp3);
    
    B0 = zeros(size(tmp1));
    [mm,nn]=size(tmp1);
    B1 = sparse(i1,j1,ones(size(i1)),mm,nn);
    B2 = sparse(i2,j2,ones(size(i2)),mm,nn);
    B3 = sparse(i3,j3,ones(size(i3)),mm,nn);
    B = B0+B1+2*B2+3*B3 + 1;
    cm = [0 0 0;1 0 0;0 1 0;0 0 1];
    figure(2);
    clf;
    hold off;
    plot(j1/15,i1,'b.','MarkerSize',2);
    hold on;
    plot(j2/15,i2,'g.','MarkerSize',2);
    plot(j3/15,i3,'y.','MarkerSize',2);
    axis off;
    xlabel('Matching vectors');
    ylabel('Channels');
    printFig('ransac');
    
    
    rowsel = (rand(1,nn)>0.6);
    rowsel(bres.pp)=ones(1,6);
    figure(2); clf;
    colormap(cm);
    hold off;
    image(B(:,find(rowsel)));
    axis off;
    %xlabel('Matching vectors');
    %ylabel('Channels');
    print -dpdf ransac.pdf
    %printFig('ransac');
    
    
end

%% Estimate r,s and refine o (PART C)

partcsavename = [savefolder expname 'C.mat'];
if ~exist(partcsavename),
    [cmatches,cres] = estimate_sr_from_matcheso(bmatches,o,bres,settings);
    eval(['save ' partcsavename ' cmatches cres']);
else
    eval(['load ' partcsavename]);
end

cres = myplot(rns,cmatches,cres,settings);


%%
if 0,
    %%
    %x = blubb.eres.x;
    x = x;
    u = rawmatches.u;
    u = u(:,:);
    [y,o,inl]=tdoa_trilateration_y(u,x);
    y = real(y);
    o = real(o);
    figure(10); clf; hold off;
    plot(sum(inl));
    [I,J]=find(inl);
    D = u(find(inl));
    [res,Dproj]=calcres(D,I,J,x,y,o);
    resm = full(sparse(I,J,res));
    figure(11);
    plot(resm','.')
    if 1,
        dres.x = x;
        dres.y = y;
        dres.o = o;
        dres.inliers = inl;
        dres.res =res; 
        dres.resm = resm;
        dres.xinorm = [1 2 4];
        dmatches = rawmatches;
        dmatches.uinliers = inl;
    end
    %%
    for kk = 2:8;
    figure(5); colormap(gray); clf;
    imagesc(settings.xtid,settings.ymeter,scores{1,kk});
    hold on;
    plot(dmatches.utimes,dmatches.u(kk,:)-dmatches.u(1,:),'b.'); 
    dproj = sqrt(sum( (y-repmat(x(:,kk),1,1218)).^2 ) ) - sqrt(sum( (y-repmat(x(:,1),1,1218)).^2 ) );
    plot(dmatches.utimes,dproj,'g.'); 
    title(num2str(kk));
    pause;
    end;
    %%
end
if 0,
    jj = 1;
   [yy,oo,inl_one_y,nr_inlierid,err_rms] = tdoa_trilateration_y_one_ransac(...
       u(:,jj),x,10,0.001,1:8);
   [I,J]=find(ones(8,1));
   [res,Dproj]=calcres(u(:,jj),I',J',x,yy,oo);

end

%% Try to find more inliers among the remaining columns in uij




partdsavename = [savefolder expname 'D.mat'];
if ~exist(partdsavename),
    x = cres.x
    u = rawmatches.u;
    [y,o,inl]=tdoa_trilateration_y(u,x);
    [I,J]=find(inl);
    D = u(find(inl));
    [res,Dproj]=calcres(D,I,J,x,y,o);
    resm = full(sparse(I,J,res));
    figure(6);
    plot(res,'.');
    ok = abs(res)<0.05;
    okm = full(sparse(I,J,ok));
    okj = (sum(okm)>=5);
    fokj = find(okj);
    % Packa in i dmatches
    dmatches.uij = rawmatches.uij(:,fokj);
    dmatches.u = rawmatches.u(:,fokj);
    dmatches.uindex = rawmatches.uindex(:,fokj);
    dmatches.utimes = rawmatches.utimes(:,fokj);
    dmatches.uok = rawmatches.uok(:,fokj);
    dmatches.uinliers = okm(:,fokj);
    inl = dmatches.uinliers;
    u = dmatches.u;
    [I,J]=find(inl);
    D = u(find(inl));
    y = y(:,fokj);
    o = o(:,fokj);
    
    % Packa in i dres
    dres.x = x;
    dres.y = y;
    dres.o = o;
    dres.xinorm = settings.xinorm;
    myplot(rns,dmatches,dres,settings);
    x = real(dres.x);
    y = real(dres.y);
    o = real(dres.o);
    [xx4,yy4,oo4,res,jac]=bundletdoa(D,I,J,x,y,o,1,dres.xinorm);
    dres.x = xx4;
    dres.y = yy4;
    dres.o = oo4;
    myplot(rns,dmatches,dres,settings);
    eval(['save ' partdsavename ' dmatches dres']);
else
    eval(['load ' partdsavename]);
end

dres = myplot(rns,dmatches,dres,settings);

%% Plot

if settings.doplot,
    figure(5); clf;
    hist(dres.res(:)*1000,100);
    %title(['Histogram of residuals']);
    xlabel('Residuals [mm]');
    ylabel('Frequency');
    printFig('residuals');
end;


if settings.doplot,
    figure(2);
    clf; hold off;
    rita3b(cres.y,'r.');
    hold on
    rita3b(cres.y,'r-');
    rita3b(cres.x,'b.');
    for mi = 1:settings.mm;
        hh = text(cres.x(1,mi)+0.05,cres.x(2,mi)+0.05,cres.x(3,mi)+0.05,num2str(mi));
        set(hh,'FontSize',30);
    end
    %xlabel('x-axis [m]');
    %ylabel('y-axis [m]');
    %zlabel('z-axis [m]');
    %printFig('3dplot');
    print -depsc 3dplot3.eps
end;

%%





%% Refinement??
%keyboard;
partesavename = [savefolder expname 'E.mat'];
if ~exist(partesavename),
    ematches = dmatches;
    for kk = 1:size(dmatches.uinliers,2);
        t0 = dmatches.utimes(kk);
        uvc = sqrt(sum( (dres.x - repmat(dres.y(:,kk),1,8)).^2 ) )' + dres.o(kk);
        uvc = uvc-uvc(1);
        % cutout;
        w = settings.wlist(1);
        
        for ii = 1:10;
            %[data0]=sampleWithU(a,t0,uvc,w,settings);
            [data,datap]=sampleWithUinterp(a,t0,uvc,w,settings);
            [u,s,v]=svd(data,0);
            %figure(1); clf;
            %plot(v(:,1)');
            dataapp = u(:,1)*s(1,1)*v(:,1)';
            err(ii)=norm(data-dataapp,'fro');
            %Flytta subpixelsteg
            dx = zeros(8,1);
            %keyboard;
            for kkk = 1:8,
                % dataapp = data + dt*datap;
                dx(kkk) = (dataapp(kkk,:)-data(kkk,:))/datap(kkk,:);
                %[f1,t1] = subpixelshift(row,ilist,t0,pattern);
            end;
            
            troskel = settings.v/settings.sr;
            troskel = troskel;
            
            if norm(dx)>troskel,
                %norm(dx)/troskel
                dx = troskel*dx/norm(dx);
            end
            
            [norm(dx) err(ii)]
            
            datapred = data + diag(dx)*datap;
            [datany,datapny]=sampleWithUinterp(a,t0,uvc+dx,w,settings);
            %norm(datapred-data)
            %norm(datapred-datany)
            %norm(data-dataapp);
            
            uvc = uvc+dx;
            dxs(:,ii)=dx;
            %figure(2); plot(data'); axis([200 400 -0.01 0.015]);
            
        end
        err
        ematches.u(:,kk)=uvc;
    end
    
    eres = dres;
    x = eres.x;
    y = eres.y;
    o = eres.o;
    ematches.uinliers = ones(size(ematches.uinliers));
    inl = ematches.uinliers;
    u = ematches.u;
    [I,J]=find(inl);
    D = u(find(inl));
    [xx5,yy5,oo5,res5,jac5]=bundletdoa(D,I,J,x,y,o,1,eres.xinorm);
    
    eres.x = xx5;
    eres.y = yy5;
    eres.o = oo5;
    eres.res = res5;
    eres.resm = sparse(I,J,res5,size(ematches.u,1),size(ematches.u,2));
    
    ematches.uij = ematches.u*settings.sr/settings.v;
    ematches.uok = ematches.uinliers;
    
    eval(['save ' partesavename ' ematches eres']);
else
    eval(['load ' partesavename]);
end

eres = myplot(rns,ematches,eres,settings);


%%

if 1,
    eres = myplot(rns,ematches,eres,settings);
    x = eres.x;
    y = eres.y;
    o = eres.o;
    ematches.uinliers = ones(size(ematches.uinliers));
    inl = ematches.uinliers;
    u = ematches.u;
    [I,J]=find(inl);
    D = u(find(inl));
    [xx5,yy5,oo5,res5,jac5]=bundletdoasmooth(D,I,J,x,y,o,1,eres.xinorm);
    fres = eres;
    fres.x = xx5;
    fres.y = yy5;
    fres.o = oo5;
    fres = myplot(rns,ematches,fres,settings);
    fmatches = ematches;
end

%% It ends here, the rest is old code

%% Test to see if one can make a bundle adjustment with smoothness constraint

if 0,
    x = eres.x;
    y = eres.y;
    o = eres.o;
    ematches.uinliers = ones(size(ematches.uinliers));
    inl = ematches.uinliers;
    u = ematches.u;
    [I,J]=find(inl);
    D = u(find(inl));
    [xx5,yy5,oo5,res5,jac5]=bundletdoasmooth(D,I,J,x,y,o,1,eres.xinorm);
    
    eres.x = xx5;
    eres.y = yy5;
    eres.o = oo5;
    eres.res = res5(1:length(I));
    eres.resm = sparse(I,J,res5(1:length(I)),size(ematches.u,1),size(ematches.u,2));
end

%% 


if 0,
    %%
    figure(10);
    clf;
    hold off;
    for kk = 1:8;
        blubb = norm(a(kk,:));
        subplot(8,1,kk); plot(a(kk,1000:1000:end));
        %axis([0 3070976/100 -0.0001 0.0001]);
        axis off
    end;
end

if 0
    litet = 0.00001;
    dx = ones(8,1)*litet;
    [data00,datap00]=sampleWithUinterp(a,t0,uvc,w,settings);
    [data11,datap11]=sampleWithUinterp(a,t0,uvc+dx,w,settings);
    blubb1 = [ (data11-data00)/litet];
    [norm(blubb1-datap00) norm(datap00) norm(blubb1)]
    
end

if 0,
    figure(1);
    imagesc(abs([data;dataapp;data-dataapp]))
    figure(2);
    plot(data')
    figure(3);
    plot(dataapp');
end


%% Make comparisons between different reconstructions

if 0,
    savefolder = '/Users/kalle/Documents/projekt/antennaResection/data/experimentEhuset_20140127/';
    explist = {'space_cure_betweenMics-', ...
        'space_cure-', ...
        'space_cure_choir-', ...
        'space_talk-' ...
        };
    
    for kkkkk = 1:3,
        expname = explist{kkkkk};
        partdsavename = [savefolder expname 'D.mat'];
        data{kkkkk} = load(partdsavename);
        
    end;
    x1 = diag([1 1 -1])*data{1}.dres.x;
    x2 = diag([-1 1 -1])*data{2}.dres.x;
    x3 = diag([1 -1 1])*data{3}.dres.x;
    xm = (x1+x2+x3)/3;
    xd1 = x1-xm;
    xd2 = x2-xm;
    xd3 = x3-xm;
    
    xs = [x1(:) x2(:) x3(:)]
    
    figure(9);
    clf; hold off;
    rita3(x1,'b*');
    hold on
    rita3(x2,'r*');
    rita3(x3,'g*');
    for mi = 1:settings.mm;
        text(x1(1,mi),x1(2,mi),x1(3,mi),num2str(mi));
    end
    
    xd = [xd1(:); xd2(:); xd3(:)];
    
    
    %standardavvikelse - cirka 9 mm
end

%%

if 0,
    %savefolder = '/Users/kalle/Documents/projekt/antennaResection/data/experimentEhuset_20140127/';
    %explist = {'space_cure_betweenMics-', ...
    %    'space_cure-', ...
    %    'space_cure_choir-', ...
    %    'space_talk-' ...
    %    };
    load xbass
    
    xinorm = settings.xinorm;
    xt = xbass2;
    xt = xt-repmat(xt(:,xinorm(1)),1,8);
    tmp = xt(:,xinorm(2:3));
    [qq,rr]=qr(tmp);
    xt = qq'*xt;
    
    m = size(xbass2,2);
    blubb = reshape(1:(3*m),3,m);
    dontmoveindex = [blubb(:,xinorm(1))'  blubb(2:3,xinorm(2))' blubb(3,xinorm(3))'];
    EE = speye(3*m+4*n);
    EE(:,dontmoveindex) = [];


%     for kkkkk = 1:3,
%         expname = explist{kkkkk};
%         partdsavename = [savefolder expname 'D.mat'];
%         data{kkkkk} = load(partdsavename);
%         
%     end;
    x1 = diag([1 1 -1])*data{1}.dres.x;
    x2 = diag([-1 1 -1])*data{2}.dres.x;
    x3 = diag([1 -1 1])*data{3}.dres.x;
    xm = (x1+x2+x3)/3;
    xd1 = x1-xm;
    xd2 = x2-xm;
    xd3 = x3-xm;
    
    xs = [x1(:) x2(:) x3(:)]
    
    figure(9);
    clf; hold off;
    rita3(x1,'b*');
    hold on
    rita3(x2,'r*');
    rita3(x3,'g*');
    for mi = 1:settings.mm;
        text(x1(1,mi),x1(2,mi),x1(3,mi),num2str(mi));
    end
    
    xd = [xd1(:); xd2(:); xd3(:)];
    
    
    %standardavvikelse - cirka 9 mm


    m = size(xt,2);
    n = size(yt,2);
    blubb = reshape(1:(3*m),3,m);
    dontmoveindex = [blubb(:,xinorm(1))'  blubb(2:3,xinorm(2))' blubb(3,xinorm(3))'];
    EE = speye(3*m+4*n);
    EE(:,dontmoveindex) = [];

end

%%
if 0,
    aa = load
    
end
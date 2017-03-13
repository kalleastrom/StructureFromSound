%% Ground Truth Generation. Partly manual.



%% set system specific paths in a separate script
% Don't do this in this script. Run the script separately or in startup.m
% so that this script works for everyone.
% set_localpaths

tdoasystemsettings.multipolpath = '/Users/kalle/Documents/projekt/github/multipol/';
tdoasystemsettings.datapath = '/Users/kalle/Documents/projekt/github/StructureFromSound/data/';

%% add paths

addpath(tdoasystemsettings.multipolpath);
addpath([tdoasystemsettings.datapath filesep 'sfsdb']);
addpath gcctracking_rex
addpath toa
addpath tdoa
addpath sfs_systems
addpath tools

%% Obtain microphones from another run

%gtloadname = ['/Users/kalle/Documents/projekt/SFS/data/sfsgt/grieg5/gt.mat'];
%load(gtloadname);
% Extract parts of the ground truth
% gt.x = result.res.x;
% gt.x_mirrors = result.speglingar;
% gt.y = result.res.y;
% gt.t = result.matches.utimes;
% gt.x = dres.x;
% gt.x_mirrors = speglingar;

% Later we are going to save the ground truth for the new data
% gtsavename = [settings.gtDir 'gt.mat'];
% if ~exist(settings.gtDir,'dir'),
%     mkdir(settings.gtDir);
% end
% eval(['save ' gtsavename ' gt']);
%

mmm = [ ...
    0 , -0.22105932 , -0.31131638 , -0.67200956 , -1.2055209 , -1.4236101 , -1.8775508 , -1.9329987 ; ...
    0 , -0.01974818 , 0.63091298 , 0.87860974 , 0.78317116 , -0.21526238 , 1.110223e-16 , 0.87813412 ; ...
    0 , -1.3555765 , 2.220446e-16 , -1.1690118 , -0.049042016 , -1.2436546 , -1.110223e-16 , -1.291896 ];

gt.x = mmm;

%% Read audio files and initial settings

[a,settings]=read_from_sfsdb(tdoasystemsettings,9);
settings.gtDir = [tdoasystemsettings.datapath 'sfsgt' filesep settings.foldername filesep];

%% Calculate GCC and top 4 peaks
[settings,scores,u] = main_gcc_and_extract_peaks(a,settings);

%% Clean up???
settings.binSize = 3;         %inlier threshold
settings.minNbrOfInliers = 3; %min number of matching equations
[uclean,ucleank] = matchingdelays2(u,settings);

%%
mm1 = max(urefk{1,2});
mm = mm1;
for i1 = 1:7;
    for i2 = (i1+1):8;
        mm = max([mm;max(urefk{i1,i2})]);
    end
    
end
figure(1);
clf;
hold off;
plot(mm,'.')
hold on;
plot(mm1,'r.');
%%
if 1,
    skalf = settings.v/settings.sr;
    for i1 = 1:7;
        for i2 = (i1+1):8,
            
            
            figure(1); clf; hold off;
            plot(settings.xtid,skalf*uclean{i1,i2}','.');
            pause;
            
        end
        
    end
    
    
end

%% My special uclean algorithm

T = 0.02;
j = 1862;
b = [];
for i1 = 1:7,
    for i2 = (i1+1):8;
        tmp = skalf*uclean{i1,i2}(:,j);
        tmp = tmp(find(isfinite(tmp)));
        b = [b; i1*ones(size(tmp)) i2*ones(size(tmp)) tmp];
        
        
        
    end
end
i1 = 7;
i2 = 8;
j = 1862;
for i3 = 1:6;
    
end
seeds_sel = [81];
seeds_i = b(seeds_sel,1:2);
seeds_d = [0 b(seeds_sel,3)];
dd = NaN*ones(1,8);
dd(seeds_i)=seeds_d;
ok_in = zeros(1,8);
ok_in(seeds_i)=1;
ok_out = zeros(1,8);
ok_out(1)=1;
ok1 = find(ok_in(b(:,1)) & ok_out(b(:,2)));
ok2 = find(ok_out(b(:,1)) & ok_in(b(:,2)));
ok12 = [ok1];

tmp = [dd(b(ok2,1))'+b(ok2,3);dd(b(ok2,2))'-b(ok2,3)];
% Find the lowest value with lots of inliers

[maxv,minc,ind]=find_cluster(tmp,T);
if maxv>=2,
    
    
    ii = 1;
    i1 = b(ii,1);
    i2 = b(ii,2);
    tmp1 = b(find(b(:,1)==i1),2:3);
    tmp2 = b(find(b(:,1)==i2),2:3);
    tmp2(:,2)=tmp2(:,2)+b(ii,3);
    for i3 = setdiff(1:8,[i1 i2]);
        dd1 = tmp1(find(tmp1(:,1)==i3),2);
        dd2 = tmp2(find(tmp2(:,1)==i3),2);
        blubb = repmat(dd1,1,length(dd2))-repmat(dd2',length(dd1),1);
        [II,JJ]=find(abs(blubb)<T);
        if length(II)>0,
            
        end
    end
    tmp1'
    tmp2'
end


%% My special ransac using multipath components

skalf = settings.v/settings.sr;
x_direct = gt.x;;
%x_all = [ [gt.x;1:8] gt.x_mirrors];
x_all = [ [gt.x;1:8]];
ynew = [];
onew = [];
jnew = [];
inlnew = [];
kk = 0;
for j = 1:size(u{1,2},2); % time index
    j
    ys = zeros(3,0);
    os = zeros(1,0);
    nr_inliers = 0;
    evals = 0;
    counter = 0;
    for kkk = 1:1000;
        refCh = randi(8);
        selfrom = setdiff(1:8,refCh);
        othCh = selfrom(randperm(7,3));
        %Choose random distance data from u (from gcc)
        d1 = skalf*u{refCh,othCh(1)}(randi(settings.nbrOfPeaks),j);
        d2 = skalf*u{refCh,othCh(2)}(randi(settings.nbrOfPeaks),j);
        d3 = skalf*u{refCh,othCh(3)}(randi(settings.nbrOfPeaks),j);
        d123 = [0 d1 d2 d3]';
        if all(isfinite(d123)),
            [ysols,osols]=tdoa_trilateration_y_one_point(d123,gt.x(:,[refCh othCh]));
            for ii = 1:2;
                counter = counter+1;
                yy = real(ysols(:,ii));
                oo = real(osols(:,ii));
                d_proj1 = tdoa_calc_u_from_xyo(x_direct,yy,oo);
                d_proj2 = tdoa_calc_u_from_xyo(x_all(1:3,:),yy,oo);
                ys(:,counter)=yy;
                os(:,counter)=oo;
                % Complicated inlier calculations??
                nrok = 0;
                cumres = 0;
                for i1 = 1:7;
                    for i2 = (i1+1):8;
                        ok1 = find(x_all(4,:)==i1);
                        ok2 = find(x_all(4,:)==i2);
                        doa_proj = -[ ...
                            repmat(d_proj1(i1),length(ok2),1)-d_proj2(ok2);
                            d_proj2(ok1(2:end))-repmat(d_proj1(i2),length(ok1)-1,1);
                            ];
                        derr = repmat(doa_proj,1,4)-repmat(u{i1,i2}(:,j)'*skalf,length(doa_proj),1);
                        nrok = nrok + sum(min(abs(derr)) < 0.1);
                        cumres = cumres + sum( (min(abs(derr)).^2).*(min(abs(derr)) < 0.03));
                    end
                end
                %hist(nr_inliers,1:100)
                %inlierid2 = find(abs(d_proj-d) < ransac_tol);
                %keyboard;
                nr_inliers(1,counter)=nrok;
                if nrok>0,
                    evals(1,counter)= sqrt(cumres)/nrok;
                else
                    evals(1,counter) = Inf;
                end
            end
        end
    end
    [maxv,maxi]=max(nr_inliers);
    if maxv>10;
        kk = kk+1;
        ynew(:,kk) = ys(:,maxi);
        onew(:,kk) = os(:,maxi);
        jnew(:,kk) = j;
        inlnew(:,kk)=maxv;
    end
end;


%% Visualize and analyze

if 1,
    figure(1); clf; hold off,
    plot(jnew,ynew','*');
    legend({'x','y','z'})
end

if 1,
    figure(1); clf,
    colormap(gray);
    imagesc(settings.xtid,settings.ymeter,scores{3,5});
end

%%
if 1,
    skalf = settings.v/settings.sr;
    for i1 = 1:7;
        for i2 = (i1+1):8,
            
            
            figure(1); clf; hold off;
            plot(settings.xtid,skalf*u{i1,i2},'.');
            hold on;
            di1i2 = toa_calc_d_from_xy(gt.x(:,[i1 i2]),ynew);
            dd = diff(di1i2);
            plot(settings.xtid(jnew),dd,'g-');
            title([num2str(i1) ' vs ' num2str(i2)]);
            pause;
            
        end
        
    end
    
    
end


%%
tmpmeas.t = jnew;
tmpmeas.y = ynew;
tmpmeas.N = size(scores{1,2},2);
A = [eye(3) eye(3);zeros(3,3) eye(3)];
sa1 = 0.00001;
sa2 = 0.001;
sA = diag([sa1 sa1 sa1 sa2 sa2 sa2]);
lA = chol(sA)';
C = [eye(3) zeros(3,3)];
sc1 = 0.04;
sC = diag([sc1 sc1 sc1]);
lC = chol(sC)';
model.A = A;
model.sA = sA;
model.lA = lA;
model.C = C;
model.sC = sC;
model.lC = lC;
[tmpres,tmpjac,tmpxe] = calcresjac2(tmpmeas,model);

%% Try some smoothing

sigma_motion = 0.002;
sigma_measurement = 0.02;
s = (1:max(jnew));
imid = ((min(jnew)+1):(max(jnew)-1))';
IJ = [imid-1 imid imid+1]';
[res,jac,z] = calcresjac(jnew,ynew,s,IJ,sigma_motion,sigma_measurement);
figure(2);
plot(z');

asol.t = s;
asol.y = z;


%% Visualize and check for inliers

if 1,
    ojoj = [];
    TTT = 0.03;
    skalf = settings.v/settings.sr;
    for i1 = 1:7;
        for i2 = (i1+1):8,
            
            
            di1i2 = toa_calc_d_from_xy(gt.x(:,[i1 i2]),asol.y);
            dd = diff(di1i2);
            
            % Check for inlier data
            
            tmp = abs( skalf*u{i1,i2}(:,asol.t) - repmat(dd,4,1) );
            [tmpmin,tmpmini]=min(tmp);            
            ok = tmpmin < TTT;
            okt = find(ok);
            blubb = skalf*u{i1,i2}(:,okt);
            blubbi = tmpmini(okt);
            blubbi2 = sub2ind(size(blubb),blubbi,1:size(blubb,2));
            okd = blubb(blubbi2);
            
            % plot
            figure(1); clf; hold off;
            plot(settings.xtid,skalf*u{i1,i2},'.');
            hold on;
            plot(settings.xtid(asol.t),dd,'g-');
            plot(settings.xtid(okt),okd,'go');
            title([num2str(i1) ' vs ' num2str(i2)]);
            pause;

            nn = length(okt);
            ojoj = [ojoj;i1*ones(nn,1) i2*ones(nn,1) okt' okd'];
        end
        
    end
    
    %data = ojoj;
end

%% New special bundle? bundles over all pairwise tdoa measurements w_i1,i2(j)
model.smotion = 0.1;
model.smeas = 1;

asol.x = gt.x;
[asolopt,res,jac]=bundletdoasmooth_allchannels(data,asol,model,debug,xinorm);
asol = asolopt;


%% Visualize new solution and check for inliers

if 1,
    ojoj = [];
    TTT = 0.03;
    skalf = settings.v/settings.sr;
    for i1 = 1:7;
        for i2 = (i1+1):8,
            
            
            di1i2 = toa_calc_d_from_xy(asol.x(:,[i1 i2]),asol.y);
            dd = diff(di1i2);
            
            % Check for inlier data
            
            tmp = abs( skalf*u{i1,i2}(:,asol.t) - repmat(dd,4,1) );
            [tmpmin,tmpmini]=min(tmp);            
            ok = tmpmin < TTT;
            okt = find(ok);
            blubb = skalf*u{i1,i2}(:,okt);
            blubbi = tmpmini(okt);
            blubbi2 = sub2ind(size(blubb),blubbi,1:size(blubb,2));
            okd = blubb(blubbi2);
            
            % plot
            figure(1); clf; hold off;
            plot(settings.xtid,skalf*u{i1,i2},'.');
            hold on;
            plot(settings.xtid(asol.t),dd,'g-');
            plot(settings.xtid(okt),okd,'go');
            title([num2str(i1) ' vs ' num2str(i2)]);
            pause;

            nn = length(okt);
            ojoj = [ojoj;i1*ones(nn,1) i2*ones(nn,1) okt' okd'];
        end
        
    end
    
    %data = ojoj;
end



%%
if 1,
    jstart =  789;
    jstop  = 2763;
    jall = jstart:jstop;
    dimi = 3;
    figure(1); clf; hold off,
    plot(jnew,ynew(dimi,:),'*');
    axis([jstart jstop -3 0]);
    [xx,yy]=ginput();
    
    tmp = interp1(xx,yy,jall)
    
    zall = tmp;
    yall = tmp;
    xall = tmp;
    man_y = [xall;yall;zall];
    
    
    skalf = settings.v/settings.sr;
    for i1 = 1:7;
        for i2 = (i1+1):8,
            
            
            figure(1); clf; hold off;
            plot(settings.xtid,skalf*u{i1,i2},'.');
            hold on;
            di1i2 = toa_calc_d_from_xy(gt.x(:,[i1 i2]),man_y);
            dd = diff(di1i2);
            plot(settings.xtid(jall),dd,'g-');
            title([num2str(i1) ' vs ' num2str(i2)]);
            pause;
            
        end
        
    end
    
    
end

%%

for i1 = 1:7;
    for i2 = (i1+1):8,
        
        
        figure(1); clf; hold off;
        plot(settings.xtid,skalf*u{i1,i2},'.');
        hold on;
        di1i2 = toa_calc_d_from_xy(gt.x(:,[i1 i2]),gt.y);
        dd = diff(di1i2);
        plot(gt.t,dd,'g-');
        title([num2str(i1) ' vs ' num2str(i2)]);
        pause;
        
    end
    
end

%% Start and stop
jstart = 789
jstop = 2763;
% jstart  = 2317
% jstop = 2317


%% Remove outliers and smooth

jwid = 50;
ransac_tol = 0.01;
ransac_k = 100;
blubb1 = NaN*ones(length(jstart:jstop),length(jnew),3);
blubb2 = NaN*ones(length(jstart:jstop),length(jnew),3);
my_pred = NaN*ones(length(jstart:jstop),3);
my_err  = NaN*ones(length(jstart:jstop),3);
jall = jstart:jstop;
for j = jstart:jstop;
    jk = find((jstart:jstop)==j)
    oki = find( (jnew > j-jwid ) & (jnew < j+ jwid));
    j_interp = jall( find( (jall > j-jwid ) & (jall < j+ jwid) ) );
    %tmp = [ynew(:,oki); onew(:,oki); jnew(:,oki)]; % Hmm. The offsets are really bad
    for dimi = 1:3;
        tmp = [ynew(dimi,oki); jnew(:,oki)]; % Hmm. The offsets are really bad
        max_nr_inl = 0;
        for kk = 1:ransac_k,
            sel_points = randperm(size(tmp,2),2);
            jj = tmp(end,:);
            j1 = tmp(end,sel_points(1));
            j2 = tmp(end,sel_points(2));
            y1 = tmp(1:(end-1),sel_points(1));
            y2 = tmp(1:(end-1),sel_points(2));
            
            predy = (y2-y1)*(jj-j1)/(j2-j1) + repmat(y1,1,size(tmp,2));
            err = sqrt(sum( (predy-tmp(1:(end-1),:)).^2 ,1));
            nr_inl = sum( (abs(err)<ransac_tol));
            if nr_inl>max_nr_inl,
                max_nr_inl=nr_inl;
                best_predy = predy;
                best_err = err;
                best_sel = sel_points;
            end
        end
        %blubb1(jk,oki,dimi)=best_predy;
        %blubb2(jk,oki,dimi)=best_err;
        %my_predall = (y2-y1)*(j_interp-j1)/(j2-j1) + repmat(y1,1,size(tmp,2));
        sel_points = best_sel;
        jj = tmp(end,:);
        j1 = tmp(end,sel_points(1));
        j2 = tmp(end,sel_points(2));
        y1 = tmp(1:(end-1),sel_points(1));
        y2 = tmp(1:(end-1),sel_points(2));
        
        my_pred(jk,dimi)=(y2-y1)*(j-j1)/(j2-j1) + repmat(y1,1,1);
        
    end;
end

%%
if 0,
    figure(1); clf; hold off;
    plot(jnew,ynew','*');
    hold on
    plot(jall,my_pred,'-');
    
end

%% More smoothing? Or fitting to direct path data?

%%

if 0,
    gt.soundpiece(1).tstart = settings.xtid(jstart);
    gt.soundpiece(1).tstop = settings.xtid(jstop);
    gt.soundpiece(1).t = settings.xtid(jall);
    gt.soundpiece(1).y = man_y;
    
end
if 0,
    gt.soundpiece(1).tstart = min(gt.t);
    gt.soundpiece(1).tstop = max(gt.t);
    gt.soundpiece(1).t = gt.t;
    gt.soundpiece(1).y = gt.y;
    
end

gt.soundpiece(1).tstart = settings.xtid(jstart);
gt.soundpiece(1).tstop = settings.xtid(jstop);
gt.soundpiece(1).t = settings.xtid(jall);
gt.soundpiece(1).y = my_pred';

%%
if isfield(gt, 't')
    % Field is there.  Remove it.
    gt = rmfield(gt, 't')
end
if isfield(gt, 'y')
    % Field is there.  Remove it.
    gt = rmfield(gt, 'y')
end

gtsavename = [settings.gtDir 'gt.mat'];
if ~exist(settings.gtDir,'dir'),
    mkdir(settings.gtDir);
end
eval(['save ' gtsavename ' gt']);


%%
% [settings,scores,rawmatches,result] = main_gcctracking_rex(a,settings);
% [settings,scores,u] = main_gcc_and_extract_peaks(a,settings);
%  rns = scores(settings.refChannel,:);
%  cres.x = gt.x;
%  cmatches = [];
%  [dmatches,dres] = find_more_inliers_uij(rawmatches,rns,cmatches,cres,settings);
% u = getdelays(scores,settings);  % Get the top k=4 tdoa measurments
% %%
% tmp = [0 settings.xtid];
% skalf = settings.v/settings.sr;

%%
% for i1 = 1:7;
%     for i2 = (i1+1):8,
%
%
%         figure(1); clf; hold off;
%         plot(tmp,skalf*u{i1,i2},'.');
%         hold on;
%         di1i2 = toa_calc_d_from_xy(dres.x(:,[i1 i2]),dres.y);
%         dd = diff(di1i2);
%         plot(dmatches.utimes,dd,'g*');
%         title([num2str(i1) ' vs ' num2str(i2)]);
%         pause;
%
%     end
%
% end

function errs = plot_and_compare_result(result);

%% 1 Compare receivers

%keyboard;
x = result.res.x;
xgt = result.settings.gt.x;
[rmse,T,x_transformed]=rigidTransform(x,xgt)

figure(1); clf;
hold off;
rita3(xgt,'go');
hold on;
rita3(x_transformed,'b*');
title(['Green - Ground Truth, Blue - Estimated. RMSE: ' num2str(rmse)]);
errs.xrmse = rmse;

%% 2 Compare Sound Path

y = result.res.y;
t = result.matches.utimes;
ygt = result.settings.gt.soundpiece(1).y;
tgt = result.settings.gt.soundpiece(1).t;
y_transformed = pflat(T*fyllmedettor(y));
y_transformed = y_transformed(1:3,:);

figure(2); clf;
hold off;
rita3(xgt,'go');
hold on;
rita3(x_transformed,'b*');
rita3(ygt,'g-');
rita3(y_transformed,'b-');
title(['Green - Ground Truth, Blue - Estimated. RMSE: ' num2str(rmse)]);

%% 3 Compare Sound Path

figure(3); clf;
strs = {'x','y','z'};
for i = 1:3;
    subplot(3,1,i);
    hold off;
    plot(tgt,ygt(i,:),'go');
    hold on
    plot(t,y_transformed(i,:),'b*');
    title([strs{i} '-coordinate' ]);
end


%% 4 Compare matches

if 1,
    u = getdelays(result.scores,result.settings);
    skalf = result.settings.v/result.settings.sr;
    xtid = result.settings.xtid;
    xtid = [0 xtid];
    ymeter = result.settings.ymeter;
    y = result.res.y;
    x = result.res.x;
    ygt = result.settings.gt.soundpiece(1).y;
    xgt = result.settings.gt.x;
    speglingar = result.speglingar;
    speglingar = [[x;1:8] result.speglingar];
    speglingargt = result.settings.gt.x_mirrors;
    speglingargt = [[xgt;1:8] speglingargt];
    t = result.matches.utimes;
    tgt = result.settings.gt.soundpiece(1).t;
    t1 = min(t);
    t2 = max(t);
    t1gt = min(tgt);
    t2gt = max(tgt);
    t1o = max(t1,t1gt);
    t2o = min(t2,t2gt);
    errs.overlap = (t2o-t1o)/(t2gt-t1gt);
    okt = find( (t>=t1o) & (t<= t2o) );
    oktg = zeros(size(okt));
    okd = zeros(size(okt));
    for k = 1:length(okt);
        [tmp,tmpi]=min( abs(tgt-t(okt(k))) );
        oktg(k)=tmpi;
    end
    dtmp = toa_calc_d_from_xy(speglingar(1:3,:),y);
    dtmpgt = toa_calc_d_from_xy(speglingargt(1:3,:),ygt);
    td   = dtmp(2:8,okt)-repmat(dtmp(1,okt),7,1);
    tdgt = dtmpgt(2:8,oktg)-repmat(dtmpgt(1,oktg),7,1);
    tderr = sqrt( sum( (td-tdgt).^2 ) );
    errs.tdok = sum( abs(tderr)<0.05 )/length(tderr);  
end





%% 5 Compare matches. Delete?

if 1,
    u = getdelays(result.scores,result.settings);
    skalf = result.settings.v/result.settings.sr;
    xtid = result.settings.xtid;
    xtid = [0 xtid];
    ymeter = result.settings.ymeter;
    y = result.res.y;
    x = result.res.x;
    ygt = result.settings.gt.soundpiece(1).y;
    xgt = result.settings.gt.x;
    speglingar = result.speglingar;
    speglingar = [[x;1:8] result.speglingar];
    speglingargt = result.settings.gt.x_mirrors;
    speglingargt = [[xgt;1:8] speglingargt];
    t = result.matches.utimes;
    tgt = result.settings.gt.soundpiece(1).t;
    t1 = min(t);
    t2 = max(t);
    t1gt = min(tgt);
    t2gt = max(tgt);
    t1o = max(t1,t1gt);
    t2o = min(t2,t2gt);
    errs.overlap = (t2o-t1o)/(t2gt-t1gt);
    okt = find( (t>=t1o) & (t<= t2o) );
    oktg = zeros(size(okt));
    for k = 1:length(okt);
        [tmp,tmpi]=min( abs(tgt-t(okt(k))) );
        oktg(k)=tmpi;
    end
    dtmp = toa_calc_d_from_xy(speglingar(1:3,:),y);
    dtmpgt = toa_calc_d_from_xy(speglingargt(1:3,:),ygt);
    for i1 = 1:7;
        for i2 = (i1+1):8;
            idi1 = find(speglingar(4,1:8)==i1);
            idi2 = find(speglingar(4,1:8)==i2);
            idi1gt = find(speglingargt(4,1:8)==i1);
            idi2gt = find(speglingargt(4,1:8)==i2);
            idi1 = find(speglingar(4,:)==i1);
            idi2 = find(speglingar(4,:)==i2);
            idi1gt = find(speglingargt(4,:)==i1);
            idi2gt = find(speglingargt(4,:)==i2);
            
            figure(2);
            clf;
            colormap(gray);
            imagesc(xtid,ymeter,max(result.scores{i1,i2},-0.1))
            %plot(xtid,skalf*u{i1,i2}','r.');
            hold on;
            for k1 = idi1;
                for k2 = idi2;
                    plot(t,dtmp(k2,:)-dtmp(k1,:),'b.');
                end
            end
            for k1 = idi1gt;
                for k2 = idi2gt;
                    plot(tgt,dtmpgt(k2,:)-dtmpgt(k1,:),'g.');
                end
            end
            title(num2str([i1 i2]));
            axis([min(xtid) max(xtid) min(ymeter) max(ymeter)]);
            %axis ij;
            %keyboard;
        end
    end
    
end






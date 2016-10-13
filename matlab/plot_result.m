function plot_result(result);

%% Plot receivers and senders

x = result.res.x;
y = result.res.y;
speglingar = result.speglingar;

%takplanid = [1     4     6     8    12    14    16    19];
%takgolvid = [2     5     7     9    13    15];
% nn = mean(spegelplan(1:3,takplanid)')';
% [uuu,sss,vvv]=svd(nn);
% uuu = diag([1 1 -1])*[0 0 1;0 -1 0;1 0 0]*uuu';
% TT = [uuu [0;0;1.65];zeros(1,3) 1];
% inv(TT)'*spegelplan
%
TT = eye(4); % One could use a change of coordinates here

figure(1); clf;
rita3(pflat(TT*fyllmedettor(x)),'g*');
hold on;
rita3(pflat(TT*fyllmedettor(speglingar(1:3,:))),'r*');
rita3(pflat(TT*fyllmedettor(y)),'g-');
%title('Valt koordinatsystem så att spegelplanen är parallella med xy-planen ungefär');

%% Plot reflective planes
%
% xx = (-0.5):0.5:4;
% yy = -3:0.5:3;
% zz = 0;
% for x0 = xx;
%     plot3([x0 x0],[min(yy) max(yy)],[zz zz],'b');
% end;
% for y0 = yy;
%     plot3([min(xx) max(xx)],[y0 y0],[zz zz],'b');
% end;
%
% xx = (-0.5):0.5:4;
% yy = -3:0.5:3;
% zz = 2.8;
% for x0 = xx;
%     plot3([x0 x0],[min(yy) max(yy)],[zz zz],'b');
% end;
% for y0 = yy;
%     plot3([min(xx) max(xx)],[y0 y0],[zz zz],'b');
% end;
%
% xx = -0.5;
% yy = -3:0.5:3;
% zz = 0:0.4:2.8;
% for z0 = zz;
%     plot3([xx xx],[min(yy) max(yy)],[z0 z0],'b');
% end;
% for y0 = yy;
%     plot3([xx xx],[y0 y0],[min(zz) max(zz)],'b');
% end;
%
% axis equal
%

%%

xtid = result.settings.xtid;
xtid = [0 xtid];
ymeter = result.settings.ymeter;
y = result.res.y;
x = result.res.x;
speglingar = result.speglingar;
speglingar = [[x;1:8] result.speglingar];
%speglingar = [[x;1:8]];
yt = result.matches.utimes;
dtmp = toa_calc_d_from_xy(speglingar(1:3,:),y);
for i1 = 1:7;
    for i2 = (i1+1):8;
        idi1 = find(speglingar(4,:)==i1);
        idi2 = find(speglingar(4,:)==i2);
        
        figure(2);
        clf;
        colormap(gray);
        imagesc(xtid,ymeter,max(result.scores{i1,i2},-0.1))
        hold on;
        for k1 = idi1;
            for k2 = idi2;
                plot(yt,dtmp(k2,:)-dtmp(k1,:),'.');
            end
        end
        title(num2str([i1 i2]));
        axis([min(xtid) max(xtid) min(ymeter) max(ymeter)]);
        axis ij;
        %pause;
    end
end

%%

xtid = result.settings.xtid;
xtid = [0 xtid];
ymeter = result.settings.ymeter;
y = result.res.y;
x = result.res.x;
speglingar = result.speglingar;
speglingar = [[x;1:8] result.speglingar];
speglingar = [[x;1:8]];
yt = result.matches.utimes;
dtmp = toa_calc_d_from_xy(speglingar(1:3,:),y);
for i1 = 1:7;
    for i2 = (i1+1):8;
        idi1 = find(speglingar(4,:)==i1);
        idi2 = find(speglingar(4,:)==i2);
        
        figure(2);
        clf;
        subplot(1,2,1);
        colormap(gray);
        imagesc(xtid,ymeter,max(result.scores{i1,i2},-0.1))
        hold on;
        subplot(1,2,2);
        hold on;
        for k1 = idi1;
            for k2 = idi2;
                plot(yt,dtmp(k2,:)-dtmp(k1,:),'-');
            end
        end
        title(num2str([i1 i2]));
        axis([min(xtid) max(xtid) min(ymeter) max(ymeter)]);
        axis ij;
        pause;
    end
end



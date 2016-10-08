function [speglingar,spegelplan]=finding_multipath_components(fres,fmatches,scores,settings);
% Try to estimate mirrored receivers. Experimental

%% load in some data
%load fresfmatchesscores
speglingar = [];

%% Finding multipath components

ind = fmatches.uindex;
u = getdelays(scores,settings);  % Get the top k=4 tdoa measurments

%%

y = fres.y;
x = fres.x;
o = fres.o;

%% Since we got the same -d1 distance for every i2
% I might as well plot them in the same plot
% This is done here.

% Go through all channels.
% The idea is to search for mirrored versions of channelMirrorSearch.
for channelMirrorSearch = 1:8;
    otherChannels = setdiff(1:8,channelMirrorSearch);
    
    % Generate data for ransac
    d1 = sqrt( sum( (y-repmat(x(:,channelMirrorSearch),1,size(y,2))).^2) );
    % Distances between sound s and microphone r_{i,1} in the paper
    figure(1);
    clf; hold on;
    skalf = settings.v/settings.sr;
    data = zeros(4,length(ind),7); % 4 = top 4 peaks,
    for kkk = otherChannels;
        %    figure(kkk);
        dd = sqrt( sum( (y-repmat(x(:,kkk),1,size(y,2))).^2) );
        %    hold off;
        data(:,:,find(otherChannels==kkk))=(skalf*u{channelMirrorSearch,kkk}(:,ind)-repmat(dd,4,1));
        plot(ind,(skalf*u{channelMirrorSearch,kkk}(:,ind)-repmat(dd,4,1))','.');
        hold on;
        plot(ind,-d1,'g-');
        %pause;
    end
    
    
    % It looks ok. There are several votes for the direct path
    % Now remove some of the noise.
    % Keep only those, where there are many votes for a distance
    % measurement.
    % Check how many votes that are close to a given measurement.
    nrinl = zeros(size(data));
    inliddata = zeros(size(data));
    TTT = 0.01; % Inlier threshold 1 cm
    for jj = 1:length(ind);
        %jj = ind(jjk);
        tmp = data(:,jj,:);
        tmp = tmp(:);
        for ii = 1:4;
            for kk = 1:7;
                nrinl(ii,jj,kk) = sum(  (abs(tmp-data(ii,jj,kk)) < TTT ));
                inliddata(ii,jj,kk)=ind(jj);
            end;
        end
    end
    
    % Also remove original track - the direct path
    for jj = 1:length(ind);
        for ii = 1:4;
            for kk = 1:7;
                if abs(data(ii,jj,kk)+d1(jj))<0.03, % Remove everything closer than 3 cm
                    nrinl(ii,jj,kk)=0;
                end
            end;
        end
    end
    
    fortsatt = 1;
    
    while fortsatt,
        %
        ok = find(nrinl>=3); % Remove those that have too few votes
        % figure(1); clf;
        % plot( inliddata(ok), data(ok),'.');
        % Looks much cleaner doesn't it
        
        % Run ransac over TOA3D
        
        % set up data for ransac
        ok = find(nrinl>=3);
        iii = inliddata(ok); % id nr
        ddd = -data(ok)';
        converter = zeros(1,max(ind));
        converter(ind)=1:length(ind);
        iii = converter(iii);
        nnn = length(ddd);
        
        %
        % Ransac
        maxnrinls = 0;
        bestx = zeros(3,1);
        bestinls = zeros(size(ddd));
        for tryid = 1:100;
            % Pick 3 random data points
            blubb = randperm(nnn);
            blubb = blubb(1:3);
            % data is
            %ddd(blubb)
            %y(:,iii(blubb))
            if length(unique(iii(blubb)))==3, % Check that the points are different
                % Calculate mirrored receiver using trilateration. -> two
                % solutions
                xsols=toa_trilateration_one_point(ddd(blubb)',y(:,iii(blubb)));
                if isreal(xsols),
                    for kkk = 1:2; % Check both solutions
                        % Calculate reprojected path
                        dproj =  sqrt( sum( ( y(:,iii)-repmat(xsols(:,kkk),1,nnn) ).^2 ) );
                        
                        
                        nrinls = sum( abs(dproj-ddd)<0.03); % Check for inliers (3 cm off)
                        if nrinls>maxnrinls,
                            nrinls
                            maxnrinls=nrinls;
                            bestx = xsols(:,kkk);
                            bestinls = find(abs(dproj-ddd)<0.03);
                        end
                        
                        if 0,
                            figure(1); clf;
                            hold off;
                            plot(iii,ddd,'b.');
                            hold on;
                            plot(iii,dproj,'g.');
                            title(num2str(nrinls));
                            %pause;
                        end
                        
                    end
                end
            end
        end;
        
        %
        dproj =  sqrt( sum( ( y(:,iii)-repmat(bestx,1,nnn) ).^2 ) );
        nrinls = sum( abs(dproj-ddd)<0.03);
        if 1,
            figure(1); clf;
            hold off;
            plot(iii,ddd,'b.');
            hold on;
            plot(iii,dproj,'g.');
            title(num2str([channelMirrorSearch nrinls]));
            pause(0.1);
            pause;
        end
        
        %
        if nrinls>100, % Make some sort of criterion for it being ok.
            speglingar = [speglingar [bestx;channelMirrorSearch]];
            dproj2 = sqrt( sum( (y-repmat(bestx(:,1),1,size(y,2))).^2) );
            
            % remove these inliers for further search
            for jj = 1:length(ind);
                for ii = 1:4;
                    for kk = 1:7;
                        if abs(data(ii,jj,kk)+dproj2(jj))<0.03,
                            nrinl(ii,jj,kk)=0;
                        end
                    end;
                end
            end
            
            
            figure(4); clf;
            rita3(x,'g*');
            hold on;
            rita3(speglingar(1:3,:),'r*');
            rita3(y,'g-');
            fortsatt = 1;
        else
            fortsatt = 0;
        end;
        
        %
        
        save speglingartmp speglingar
    end
    
end;


%% Hitta spegelplan
spegelplan = speglingar;
for ii = 1:size(speglingar,2);
    pid = speglingar(4,ii);
    p1 = x(:,pid);
    p2 = speglingar(1:3,ii);
    p12 = (p1+p2)/2;
    n = p2-p1;
    n = n/norm(n);
    plan = [n;-n'*p12];
    spegelplan(1:4,ii)=plan;
end


%% Cluster the reflective planes
%[aaa,bbb]=kmeans(spegelplan',6);
%spegelplan(:,find(aaa==1))

% takplanid = [1     4     6     8    12    14    16    19]; %Manual identification
% takgolvid = [2     5     7     9    13    15];  %Manual identification
% nn = mean(spegelplan(1:3,takplanid)')';
% [uuu,sss,vvv]=svd(nn);
% uuu = diag([1 1 -1])*[0 0 1;0 -1 0;1 0 0]*uuu';
% TT = [uuu [0;0;1.65];zeros(1,3) 1];
% inv(TT)'*spegelplan
%
% figure(5); clf;
% rita3(pflat(TT*fyllmedettor(x)),'g*');
% hold on;
% rita3(pflat(TT*fyllmedettor(speglingar(1:3,:))),'r*');
% rita3(pflat(TT*fyllmedettor(y)),'g-');

%%




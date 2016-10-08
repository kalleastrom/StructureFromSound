function [o,v,inliers,nr_inliers,err]=tdoa_offset_one_ransac(utmp,u_anchor,o_anchor,Abase,rankk,ransac_k,ransac_tol,inlierid);
% [y,inliers,nr_inliers,err_rms]=tdoa_offset_one_ransac(d,x,ransac_k,ransac_tol,inlierid);
% trilateration of one point using ransac
%
%

if nargin<3,
    ransac_k = 10;
end

if nargin<4,
    ransac_tol = Inf;
end

if nargin<5,
    inlierid = find(isfinite(u));
end

mm0 = size(utmp,1);

utmp = utmp(inlierid);
u_anchor = u_anchor(inlierid,:);
Abase = Abase(inlierid(2:end)-1,:);

[mm,nn]=size(Abase);
%md = size(d,1);

if (rankk)==mm,
    error('Underconstrained problem - to little data to perform ransac');
else
    % Main part (more data (mm) than dimensionality of problem (rankk +1)
    
    counter = 0;
    os = zeros(1,ransac_k);
    evals = zeros(1,ransac_k);
    nr_inliers = zeros(1,ransac_k);
    vtmps = zeros(rankk,ransac_k);
    inlierlist = zeros(mm0,ransac_k);
    
    for kk = 1:ransac_k
        mysel = [1 randperm(mm-1)+1];
        mysel = mysel(1:(rankk+2));
        u_anchor0 = u_anchor(mysel);
        Abase0 = Abase(mysel(2:end)-1,:);
        utmp0 = utmp(mysel);
        [otmp0,vtmp0]=tdoa_offset_one(utmp0,u_anchor0,Abase0);
        % Reconstruct u from solution
        a1 = otmp0^2-o_anchor^2;
        br = Abase*vtmp0;
        a2end = br+a1;
        ar = [a1;a2end];
        utmpr = sqrt(ar + u_anchor.^2)+otmp0;
        %utmpr - utmp
        if 0,
            %258   168    24   106    79   158
            ll = 2;
            jj = pp(ll);
            utmp = u(:,jj);
            mysel = [1 2 3 5 6];
            u_anchor0 = u_anchor(mysel);
            Abase0 = Abase(mysel(2:end)-1,:);
            utmp0 = utmp(mysel);
            [otmp0,vtmp0]=tdoa_offset_one(utmp0,u_anchor0,Abase0);
            otmp0
            oopt2(ll)
            
            % Reconstruct data from Abase0, vtmp0
            Abase*vtmp0
            a=(utmp-otmp0).^2 - u_anchor.^2
            b=a(2:end)-a(1)
            % Now Abase*vtmp0 = b at least at positions
            % mysel(2:end)
            
            % a(1) = otmp0^2-oopt2(1)^2;
            a1 = otmp0^2-oopt2(1)^2;
            br = Abase*vtmp0;
            a2end = br+a1
            ar = [a1;a2end];
            utmpr = sqrt(ar + u_anchor.^2)+otmp0;
            utmpr - utmp
            
        end
        
        counter = counter+1;
        os(:,counter)=otmp0;
        vtmps(:,counter)=vtmp0;
        inlierid0 = find(abs(utmpr-utmp) < ransac_tol);
        inlierlist(inlierid(inlierid0),counter)=ones(length(inlierid0),1);
        nr_inliers(1,counter)=length(inlierid0);
        evals(1,counter)= sqrt(sum( (utmpr(inlierid0)-utmp(inlierid0)).^2 ));
    end
    
    % Select the solution with (1) the most inliers and
    % among these the one with the lowest error score.
    [maxv,maxi]=max(nr_inliers);
    ok1 = find(nr_inliers==maxv);
    [minv,mini]=min(evals(ok1));
    besti = ok1(mini);
    
    o = os(:,besti);
    v = vtmps(:,besti);
    err = evals(1,besti);
    
    inliers = inlierlist(:,besti);
    nr_inliers=sum(inliers);
    
end


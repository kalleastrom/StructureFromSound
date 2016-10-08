function [o,v,u,A,inliers_ut,pp,cc]=tdoa_offset_ransac(u,ok);
% [o,inliers,u]=tdoa_offset_ransac(u);
% RANSAC based solver for finding offsets o
% such that double compaction of u with offsets o
% has close to rank 2
%    u(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 )) - o(j)
%
% In:
%    u       - mxn matrix with (m>=7) distances
%    ok      - logical matrix that indicates which elements are to be
%              considered. If this is not supplied, it is assumed that
%              all elements that are not NaN are ok.
%
% Out:
%    o       - 1xn matrix
%    inliers - logial matrix with inliers.
%
% Idea
%  * select random 6 columns
%  * Solve minimal case (7x6)
%  * Extend to the remaining microphones (mx6)
%  * Test to see how many of the remaining points can be reconstructed
%    within certain accuracy

[m,n]=size(u);
if nargin<2,
    ok = ~isnan(u);
end;

load options76;
ransac_k1 = 5; % How many iterations in the RANSAC loop
ransac_k2 = 5;
ransac_tol = 0.01;
rank = 3;

maxnrinliers = 0;
for kk = 1:ransac_k1
    okm = find(sum(ok)==m);
    pp = randperm(length(okm));
    pp = okm(pp(1:6));
    cc = [1 randperm(m-1)+1];
    cc = cc(1:7);
    %keyboard;
    [sols,minerr] = tdoa_offset(u(cc,pp),options76);
    if size(sols,2)>0,
        nrsols = size(sols,2);
        if nrsols >= 1,
            for ii = 1:nrsols,
                otmp = sols(:,ii);
                if isreal(otmp),
                    if all(otmp<0),
                        if 0,
                            us = u(cc,pp);
                            d = us-repmat(otmp',7,1);
                            d2 = d.^2;
                            cl = [-ones(6,1) eye(6)];
                            cr = [-ones(5,1) eye(5)];
                            svd(cl*d2*cr')
                            
                        end
                        
                        %keyboard;
                        [oopt,Aopt,paramopt,resopt1,jacopt]=bundle_rank1(u(cc,pp),ones(size(u(cc,pp))),otmp,3);
                        % Extend to the remaining row(rows)
                        [oopt2,Aopt,paramopt,resopt2,jacopt]=bundle_rank1(u(:,pp),ones(size(u(:,pp))),oopt,3);
                        inliers_d = zeros(size(ok));
                        
                        % Must remember anchor points
                        u_anchor = u(:,pp(1))-oopt2(1);
                        o_anchor = oopt2(1);
                        [uuu,sss,vvv]=svd(Aopt);
                        uuu = uuu(:,1:3);
                        vvv = sss(1:3,1:3)*vvv(:,1:3)';
                        Abase = uuu;
                        
                        
                        %inliers_d(:,pp)=ones(m,length(pp));
                        nr_inliers = 0;
                        o2 = NaN*ones(1,n);
                        v2 = NaN*ones(3,n);
                        for jj = 1:n,
                            [o20,v20,inliers0,nr_inliers0,err0]= ...
                                tdoa_offset_one_ransac(u(:,jj),u_anchor, ...
                                o_anchor,Abase,3,ransac_k2,ransac_tol,find(ok(:,jj)));
                            if nr_inliers0>=6,
                                nr_inliers = nr_inliers+1;
                                o2(:,jj)=o20;
                                v2(:,jj)=v20;
                                inliers_d(:,jj)=inliers0;
                            else
                                inliers_d(:,jj)=zeros(m,1);
                            end;
                        end
                        if nr_inliers>maxnrinliers,
                            maxnrinliers = nr_inliers;
                            %besto = otmp;
                            besto = o2;
                            bestv = v2;
                            bestu = Abase;
                            besta = Abase*v2;
                            bestpp = pp;
                            bestcc = cc;
                            bestinliers = inliers_d;
                        end
                    end
                end
            end
        end
    end
end
if maxnrinliers == 0,
    error('ojojoj');
    o = NaN;
    v = NaN;
    u = NaN;
    A = NaN;
    pp = NaN;
    cc = NaN;
    inliers_ut = bestinliers;
else
    o = besto;
    v = bestv;
    u = bestu;
    A = besta;
    pp = bestpp;
    cc = bestcc;
    inliers_ut = bestinliers;
end;
function [ematches,eres] = sub_pixel_refinement(a,dmatches,dres,settings);
% [ematches,eres] = find_more_inliers_uij(a,dmatches,dres,settings);


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

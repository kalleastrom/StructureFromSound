function slask=rita3(bild,st)
if nargin == 1,
    st='-';
end;
if 0,
    t=linspace(0,1,size(bild,2));
    %clf;
    for j=1:size(bild,2);
        plot3(bild(1,j),bild(2,j),bild(3,j),'Color',[t(j) 0 1-t(j)],'MarkerSize',100,'Marker','.');hold on;
    end;
end;
%hold off;
if 1,
    plot3(bild(1,:),bild(2,:),bild(3,:),st);
    slask=[];
end;
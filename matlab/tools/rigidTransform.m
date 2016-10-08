function [rmse T,pHatTransformed]=rigidTransform(pHat,p)
%RIGIDTRANSFORM translates and rotates points p  to align with pHat
%
%[T RMS]=rotSolOpt(p,pHat) rotates and translates the point cloud p to 
%align with pHat, omptimally in least squares sense, thus minimizing the
%RMS error SUM_i (sum_i ( ||pHat_i - (Rp_i + t) ||^2 ), where R is a rotation
%(and/or mirroring for this alg seems to work) and t is a translation.
%
%%Waning: Optimality in least squares sense only guaranteed in 3-D.
%Allthough the alg works in any dimension, I need to go through Arun and
%see that it extend to general dimensions. At first glance it looks
%optimal, but needs to be verified.
%
%Based on Arun et al., "Least-Squares Fitting of Two 3-D Point Sets"
%Code by Siomn Burgess, 2012-08-24
%
%INPUT
%   pHat    d x m   Matrix with m d-D coordinates
%   p       d x m   Matrix with m d-D coordinates
%
%OUTPUT
%   T               d+1 x d+1   Transformation matrix to be applied to pHat points (in
%                               homogenous coordinates) to align pHat with p
%   rmse                        scalar  the RMS of error, i.e. sqrt(  (SUM_i ( ||pHat_i - (Rp_i + t)||^2)/nbrPoints  )
%   pHatTransformed  d x m      Matrix with m d-D coordinates, pHat aligned
%                               with p


%%
m=size(pHat,2);
m2=size(p,2);
dim=size(pHat,1);

if m ~=m2 
    disp('not same number of points in point clouds')
    keyboard
end

% if size(pHat,2) ~= 3 || size(p,2) ~=3
%       disp('not 3 points in 3-D')
%     keyboard
% end

%centered points q
pHatMean=mean(pHat,2);
pMean=mean(p,2);
qHat = pHat-repmat(pHatMean,1,m);
q = p-repmat(pMean,1,m);

H= qHat*q';

[U,S,V] = svd(H);

R = V*U';

rmse=sqrt(sum(sum((R*qHat-q).^2))/m);

%% to make T
transl_pHat = eye(dim+1,dim+1);
transl_p = eye(dim+1,dim+1);
R_new = eye(dim+1,dim+1);

transl_pHat(1:dim, dim+1) = -pHatMean; %transpose not needed. incy wincy bit faster?
R_new(1:dim, 1:dim) = R;
transl_p(1:dim, dim+1) = pMean;

T = transl_p * R_new * transl_pHat; 

pHatTransformed = T*[pHat; ones(1,m)];
pHatTransformed = pHatTransformed(1:end-1,:);

end
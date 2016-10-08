function [y,alpha]=psphere(x);
% [y,alpha]=psphere(x) - normalization of projective points.
% Do a normalization on the projective
% points in x. Each column is considered to
% be a point in homogeneous coordinates.
% INPUT :
%  x     - matrix in which each column is a point.
% OUTPUT :
%  y     - result after normalization.
%  alpha - depth

[a,n]=size(x);
alpha=sqrt(sum(x.^2));
y=x./(ones(a,1)*alpha);

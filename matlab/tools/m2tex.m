function s = m2tex(M,npre);

if nargin<2,
npre = 2;
end
[m,n]=size(M);
disp('\begin{pmatrix}');
for i=1:m;
 row = '';
 for j=1:n-1,
  row = [row num2str(M(i,j),npre) ' & '];
 end
 row = [row num2str(M(i,n),npre) ' \\'];
 disp(row);
end

disp('\end{pmatrix}');

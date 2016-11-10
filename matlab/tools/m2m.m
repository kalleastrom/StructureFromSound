function s = m2m(M);


[m,n]=size(M);
disp('M = [ ...');
for i=1:m;
 row = '';
 for j=1:n-1,
  row = [row num2str(M(i,j),8) ' , '];
 end
 if i<m,
  row = [row num2str(M(i,n),8) ' ; ...'];
 else
  row = [row num2str(M(i,n),8) ' ];'];
 end;
 disp(row);
end

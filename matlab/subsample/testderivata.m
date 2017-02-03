% yinvers

B = [1 0 0;0 0 0; 0 0 0];
Y0 = rand(3,3);
litet = 0.0001;

Y1 = Y0 + B*litet;

% approximativ derivata
dy1 = (inv(Y1)-inv(Y0))/litet
dy2 = -inv(Y0)*B*inv(Y0)

%%
A = rand(3,3);
b = rand(3,1);

EE = eye(9);
for k = 1:9,
    B{k} = reshape(EE(:,k),3,3);
    B{k}
end

for k = 1:9,
    dinvAdx{k} = -inv(A)*B{k}*inv(A);
    dinvAdx{k}
end

tmp = zeros(3,9);
for k = 1:9,
    dinvAbdx{k} = -inv(A)*B{k}*inv(A)*b;
    dinvAbdx{k}
    tmp(:,k)=dinvAbdx{k};
end


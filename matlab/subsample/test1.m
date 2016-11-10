%% Test of subsample methods

x = 0:1000;
y = F0(x,2);
figure(1); clf;
plot(x,y);

%% Test if scale parameter works

x = -100:1:2100;

aa = 20;
y0 = F0(x,0);
ys = F0(x,aa);
xxx = -(10*aa):1:(10*aa);
ggg = 1/sqrt(2*pi*aa^2)*exp( - xxx.^2/(2*aa^2) );
ysa = conv(y0,ggg,'same');

norm(y0-ysa)
norm(ys(500:1500)-ysa(500:1500))

figure(1);
clf;
plot(xxx,ggg);

figure(2); 
clf;
plot(x,y0,'g');
hold on
plot(x,ys,'r:');
plot(x,ysa,'b--');

%% 

tau = 3.7;
x = 1:1000;
f0 = F0(x);
f1 = F0(x+tau);

%% Discrete search

tt = -10:10;
xmid = 100:900;
for k = 1:length(tt);
    err(k) = norm(f0(xmid) - f1(xmid-tt(k)));
end

[minv,mini]=min(err);
tau0 = tt(mini);

%% Test interpolation

tao0 = 3.7;
f1t = zeros(size(f1));
a = 2;
for k = xmid;
    f1t(k) = interp1d(f1,k-tau,a);
end
for k = xmid;
    f0t(k) = interp1d(f0,k,a);
end

errny1 = norm(f0(xmid)-f1t(xmid))
errny2 = norm(f0t(xmid)-f1t(xmid))

%%

%[f,dfdx]=interp1d_with_derivative(f0,x0,a);

tao0 = 3.6;
f1t = zeros(size(f1));
f1td = zeros(size(f1));
a = 2;
for k = xmid;
    [f1t(k),f1td(k)] = interp1d_with_derivative(f1,k-tau0,a);
end
for k = xmid;
    f0t(k) = interp1d(f0,k,a);
end

[f0t(xmid);f1t(xmid);f1td(xmid)]

[(f0t(xmid)-f1t(xmid))./f1td(xmid)]
%delta_tau = ???

%tau0 = tao0 + delta_tau;
% iterera här för att minimera felet med avseende på tau0

errny1 = norm(f0(xmid)-f1t(xmid))
errny2 = norm(f0t(xmid)-f1t(xmid))


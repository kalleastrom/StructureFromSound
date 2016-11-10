function y = F0(x,s);
% 

if nargin<2,
    s = 0;
end

if 0,
    n1 = 5;
    m1 = 1000*rand(1,n1);
    s1 = (10+10*rand(1,n1));
    a1 = 10*rand(1,n1);
    n2 = 5;
    a2 = 10*rand(1,n2)
    f2 = 0.1*rand(1,n2);
    th2 = (2*pi*rand(1,n2));
    save ettf0 n1 m1 s1 a1 n2 a2 f2 th2
else
    load ettf0
    s1 = sqrt(s1.^2 + s^2);
    %a1 =
    %m1 = [852.97445 , 811.35829 , 924.7631 , 133.21626 , 622.77312 ];
    %s1 = [19.490911 , 17.540713 , 12.358983 , 14.166085 , 16.034074 ];
end;

y = 0;
for k = 1:n1;
    y = y + a1(k)*(1/sqrt(2*pi*s1(k)^2))*exp( -(x-m1(k)).^2 / (2*s1(k).^2 ) );
end;

% exp(-t^2) -> sqrt(pi)*exp(-w^2/4)
% scale with a = 1/(sqrt(2)*s)
% exp(-(a*t)^2 ) = exp( -t^2/(2*s^2) -> sqrt(2)*s*sqrt(pi)*exp(-(2*s^2*w^2)/4 )
% --> (1/sqrt(2*pi*s^2)) * sqrt(2)*sqrt(pi)*s*exp(-(2*s^2*w^2)/4 )
% --> exp(-(2*s^2*w^2)/4 )

h2 = exp(-(2*s^2*f2.^2)/4 );
 for k = 1:n1;
     y = y + h2(k)*sin(f2(k)*x+th2(k));
 end;

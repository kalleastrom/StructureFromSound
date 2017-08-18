%% test for real data where translation, doppler and amplitude is estimated
% The function find_translation_doppler_amplitude does not yet use
% convolution for the implementation.


% addpath /home/gabbi/github/StructureFromSound/matlab/gcctracking_rex/
addpath ../gcctracking_rex/

%load ../../data/savefiles/tmp3/space_cure-D.mat
%skalf = settings.v/settings.sr;
%clear a
%load allt20170217

if ~exist('a'),
    % First time load in sound data a and settings, but these are quite
    % big files
    %load spacecurea; % load a real sound experiment
    %load spacecuresettings; % load settings
    load asol2
    
    u = getdelays(scores,settings);
    
end

%%
% % Study channels 1 and 2 only
% figure(1);
% plot(u{1,2}','.')
%
% % Select a shorter timeframe where the motion is
% % almost linear (with little acceleration)
% clear usel usub
% xsel = 2241:2270;
% usel = u{1,2}(xsel);
%
% figure(2);
% plot(xsel,usel);
% xt = xtid(xsel);

%%
% make a fit to the data using a low order polynomial
% Jag provade lite olika gradtal och tittade p� residualerna.
% Vid grad 3 tyckte jag att det b�rjade se bra ut.

% [p,S,mu] = polyfit(xt,usel,3);
% ucalc = polyval(p,xt,S,mu);
%
% figure(3);
% hold off;
% plot(xt,usel,'r.');
% hold on
% plot(xt,ucalc,'g.');
%
% figure(4);
% plot(xt,usel-ucalc,'.');

%%

% Convert from time to sample point and from
% meter to offset in samplepoint
% usel is already in sample points

it = 10:2843;
nup = 0000;
ndown = 1000;
resultat = zeros(3,length(it));
noise_stds = zeros(1,length(it));
skalf = settings.v/settings.sr;
mid_points = settings.dx*(it);

for kk = 1:length(it);
    kk
    ki = it(kk);
    % Calculate approximate shift using asol
    dcalc = norm(asol.y(:,ki)-asol.x(:,2)) - norm(asol.y(:,ki)-asol.x(:,1));
    ucalc = dcalc/skalf;
    % Perturb it a bit
    ucalc = round(ucalc+randn*2);
   
    %ucalc = -62
    ii = mid_points(kk);
    temp_trans = ucalc; % Use the estimated time-difference from gcc as initialization
    cutout1 = a(1,((ii-ndown):(ii+nup)));
    cutout2 = a(2,((ii-ndown):(ii+nup))+temp_trans);
    cutout1 = cutout1-mean(cutout1);
    cutout2 = cutout2-mean(cutout2);
    
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    tt = [-15 15]; % the translations to be tried
    a2 = 2;
    a22 = a2*sqrt(2);
    [z, curr_err,f0t,f1t] = find_translation_doppler_amplitude(cutout1, cutout2, thresh, a2, tt);
    z
    
    ucalc
    z(1)+temp_trans
    
    % Noise estimation
    noise_var_estimate = var(f1t-f0t) / (2*(1/sqrt(2*pi*a22^2)));
    noise_std_estimate = sqrt(noise_var_estimate);
    %
    xx = -15:15;
    %     gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
    %     tmp = conv2((f0t+f1t)/2,gg,'same');
    %     EA = 2*sum(tmp.^2);
    % re-estreck = 2*noise_std^2 * diskret deltafunction.
    % Re-esteck = 2*noise_std^2 * normalf�rdelad med std a*sqrt(2) ????
    %     geestreck = 2*noise_std_estimate^2 * (1/sqrt(2*pi*a22^2)).*exp( - (xx.^2)/(2*a22^2) );
    %     gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
    %     tmp = conv2((f0t+f1t)/2,gg,'same');
    %     tmp2 = conv2(tmp,geestreck,'same');
    %     tmp = tmp(100+1:end-100);
    %     tmp2 = tmp2(100+1:end-100);
    %     Vb = 4*sum(tmp.*tmp2);
    %     VX = Vb/(EA^2);
    %     SX = sqrt(VX)
    %     Vb_approx = noise_std_estimate^2*4*sum(tmp.*tmp);
    %     VX_approx = noise_std_estimate^2 / (sum(tmp.*tmp));
    
    if 1,
        figure(31);
        clf; hold off
        plot(f0t,'r');
        hold on;
        plot(f1t,'g');
        z
        pause(0.1);
    end
    zfix = z;
    zfix(1)=zfix(1)+temp_trans;
    noise_stds(kk)=noise_std_estimate;
    resultat(:,kk)=zfix;
    %      usub(kk)=-(curr_trans(end)-temp_trans);
    %      ustdnoise(kk)=noise_std_estimate;
    %      ustdtrans(kk)=SX;
end

%resultat(1,1:122) = resultat(1,1:122) + round(dmatches.uij(2,(it(1:122))));

%%
di1i2 = toa_calc_d_from_xy(asol.x(:,[1 2]),asol.y);
dd = diff(di1i2);
dcalc = dd(it);
ucalc = dcalc/skalf;
u_gcc = ucalc;
sj = asol.y(:,it);
r1 = asol.x(:,1);
r2 = asol.x(:,2);
nn = size(sj,2);
u_calc = sqrt(sum( (sj-repmat(r2,1,nn)).^2 ) ) - sqrt(sum( (sj-repmat(r1,1,nn)).^2 ));
d2 = sqrt(sum( (sj-repmat(r2,1,nn)).^2 ));
d1 = sqrt(sum( (sj-repmat(r1,1,nn)).^2 ));
d2d1 = d2./d1;

utimes = mid_points/settings.sr;
sjd = diff(sj')';
td = diff(utimes);
sjp = sjd./repmat(td,3,1);
n2 = (sj-repmat(r2,1,nn));
n2n = sqrt(sum( n2.^2 ));
n2 = n2./repmat(n2n,3,1);
n1 = (sj-repmat(r1,1,nn));
n1n = sqrt(sum( n1.^2 ));
n1 = n1./repmat(n1n,3,1);
n1 = n1(:,1:(end-1));
n2 = n2(:,1:(end-1));

dopp = sum(n1.*sjp)-sum(n2.*sjp);

datasel = find(data(:,1)==1 & data(:,2)==2);
ugcc_t = data(datasel,3)*1000/96000;
ugcc_d = data(datasel,4);
ugcc_u = data(datasel,4)/skalf;



%%
figure(1);
clf;
plot(ugcc_t,ugcc_d,'r.');
hold on
plot(utimes,dcalc,'g.');
plot(utimes,skalf*resultat(1,:),'b.');
xlabel('Time in seconds');
ylabel('Distance difference d2-d1 in meters');


%
figure(2); clf
plot(1.3*d2d1.^2,resultat(3,:),'.')
hold on;
plot([0 2],[0 2],'g');
%axis([0 2 0 2]);
xlabel('Distance quotient d2/d1');
ylabel('Amplitude factor');

figure(4); clf;
plot(resultat(3,:),'b.');
hold on
plot(1.3*d2d1.^2,'g');


%
figure(3); clf
plot(dopp*skalf,1-resultat(2,1:(end-1)),'.')
hold on
plot([-0.005 0.005],[-0.005 0.005])
%axis([-0.005 0.005 -0.005 0.005])

figure(6); clf;
plot(1-resultat(2,:),'b.');
hold on
plot(dopp*skalf,'g');

%%
% 
% figure(4);
% imagesc
% 
% 
% %%
% 
% [m,n]=size(im); 
% for i =1:m,
%     str = [];
%     for j = 1:n,
%         str = [str num2str(im(i,j))];
%     end
%         disp(str);
% end


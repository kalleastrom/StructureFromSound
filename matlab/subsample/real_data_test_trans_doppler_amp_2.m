%% test for real data where translation, doppler and amplitude is estimated
% The function find_translation_doppler_amplitude does not yet use
% convolution for the implementation.


% addpath /home/gabbi/github/StructureFromSound/matlab/gcctracking_rex/
addpath ../gcctracking_rex/

if ~exist('a'),
    % First time load in sound data a and settings, but these are quite
    % big files
    %load spacecurea; % load a real sound experiment
    %load spacecuresettings; % load settings
    load blubba;
    load blubbsettings;
    xtid = settings.xtid;
    ymeter = settings.ymeter;
    % Calcualte GCC-PHAT scores
    settings.v = 340;                %speed of sound
    settings.mm = size(a,1);         %number of microphones
    settings.channels = 1:2;         %channels to read
    settings.refChannel = 1;         %reference channel
    settings.nbrOfSamples = length(a);    
    %% Correlation: GCC-PHAT
    settings.wf = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function
    settings.firstSamplePoint = 1; %center sample point of first frame
    settings.frameSize = 2048;     %width of frame in sample points
    settings.dx = 1000;            %distance between frames in sample points
    settings.frameOverlap = settings.frameSize-settings.dx; %overlap between frames
    settings.sw = 800;             %clipping of search window
    scores = gccscores(a,settings);
    
    %% Delays: Find highest peaks
    settings.nbrOfPeaks = 1;       %max number of peaks
    settings.minPeakHeight = 0.01; %min value of local maxima
    %Default: [4,0.01]
    
    u = getdelays(scores,settings);
    
end

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

it = find(dmatches.uindex>2240); % Don't look a the first part, where there 
% the tracks are interupted. 
mid_points = [0 settings.isel];
nup = 0000;
ndown = 1000;
resultat = zeros(3,length(it));

for kk = 1:length(it);
kk
    ki = it(kk);
    ii = mid_points(dmatches.uindex(ki));
    temp_trans = round(dmatches.uij(2,ki)); % Use the estimated time-difference from gcc as initialization
    cutout1 = a(1,((ii-ndown):(ii+nup)));
    cutout2 = a(2,((ii-ndown):(ii+nup))+temp_trans);
    
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    tt = [-15 15]; % the translations to be tried
    a2 = 2;
    a22 = a2*sqrt(2);
    [z, curr_err,f0t,f1t] = find_translation_doppler_amplitude(cutout1, cutout2, thresh, a2, tt);
    
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
    
    resultat(:,kk)=zfix;
%      usub(kk)=-(curr_trans(end)-temp_trans);
%      ustdnoise(kk)=noise_std_estimate;
%      ustdtrans(kk)=SX;
end

%resultat(1,1:122) = resultat(1,1:122) + round(dmatches.uij(2,(it(1:122)))); 

%%
u_gcc = skalf*dmatches.uij(2,it);
sj = dres.y(:,it);
oj = dres.o(:,it);
r1 = dres.x(:,1);
r2 = dres.x(:,2);
nn = size(sj,2);
u_calc = sqrt(sum( (sj-repmat(r2,1,nn)).^2 ) ) - sqrt(sum( (sj-repmat(r1,1,nn)).^2 ));
d2 = sqrt(sum( (sj-repmat(r2,1,nn)).^2 ));
d1 = sqrt(sum( (sj-repmat(r1,1,nn)).^2 ));
d2d1 = d2./d1;

sjd = diff(sj')';
td = diff(dmatches.utimes(it));
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

%%
figure;
clf;
plot(dmatches.utimes(it),u_gcc,'r.');
hold on
plot(dmatches.utimes(it),u_calc,'g.');
plot(dmatches.utimes(it(1:384)),skalf*resultat(1,1:384),'b.');

%
figure(2); clf
plot(d2d1(1:384),resultat(3,1:384).^2,'.')

%%
figure(3); clf
plot(dopp(1:384)*skalf,1-resultat(2,1:384),'.')
hold on
plot([0 0.002],[0 0.002])
axis([0 0.002 0 0.002])


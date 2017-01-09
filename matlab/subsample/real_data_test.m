%%

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

% Study channels 1 and 2 only
figure(1);
plot(u{1,2}','.')

% Select a shorter timeframe where the motion is
% almost linear (with little acceleration)
clear usel usub
xsel = 2241:2270;
usel = u{1,2}(xsel);

figure(2);
plot(xsel,usel);
xt = xtid(xsel);

%%
% make a fit to the data using a low order polynomial
% Jag provade lite olika gradtal och tittade på residualerna.
% Vid grad 3 tyckte jag att det började se bra ut.

[p,S,mu] = polyfit(xt,usel,3);
ucalc = polyval(p,xt,S,mu);

figure(3);
hold off;
plot(xt,usel,'r.');
hold on
plot(xt,ucalc,'g.');

figure(4);
plot(xt,usel-ucalc,'.');

%%

% Convert from time to sample point and from
% meter to offset in samplepoint
% usel is already in sample points

it = settings.isel(xsel);
jt = usel;

nup = 0000;
ndown = 1000;

for kk = 1:length(xsel);
    kk
    temp_trans = usel(kk); % Use the estimated time-difference from gcc as initialization
    cutout1 = a(1,((it(kk)-ndown):(it(kk)+nup)));
    cutout2 = a(2,((it(kk)-ndown):(it(kk)+nup))+temp_trans);
    
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    tt = [-15 15]; % the translations to be tried
    a2 = 2;
    a22 = a2*sqrt(2);
    [curr_trans, curr_err,f0t,f1t] = find_translation2(cutout1, cutout2, thresh, a2, tt);
    
    % Noise estimation
    noise_var_estimate = var(f1t-f0t) / (2*(1/sqrt(2*pi*a22^2)));
    noise_std_estimate = sqrt(noise_var_estimate);
    %
    xx = -15:15;
    gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
    tmp = conv2((f0t+f1t)/2,gg,'same');
    EA = 2*sum(tmp.^2);
    % re-estreck = 2*noise_std^2 * diskret deltafunction.
    % Re-esteck = 2*noise_std^2 * normalfördelad med std a*sqrt(2) ????
    geestreck = 2*noise_std_estimate^2 * (1/sqrt(2*pi*a22^2)).*exp( - (xx.^2)/(2*a22^2) );
    gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
    tmp = conv2((f0t+f1t)/2,gg,'same');
    tmp2 = conv2(tmp,geestreck,'same');
    tmp = tmp(100+1:end-100);
    tmp2 = tmp2(100+1:end-100);
    Vb = 4*sum(tmp.*tmp2);
    VX = Vb/(EA^2);
    SX = sqrt(VX)
    Vb_approx = noise_std_estimate^2*4*sum(tmp.*tmp);
    VX_approx = noise_std_estimate^2 / (sum(tmp.*tmp));
    
    if 1,
        figure(11);
        clf; hold off
        plot(f0t,'r');
        hold on;
        plot(f1t,'g');
        pause;
    end
    
    usub(kk)=-(curr_trans(end)-temp_trans);
    ustdnoise(kk)=noise_std_estimate;
    ustdtrans(kk)=SX;
end

%%

[p,S,mu] = polyfit(xt,usub,3);
usubcalc = polyval(p,xt,S,mu);

figure(5);
hold off;
plot(xt,usel,'ro');
hold on
plot(xt,ucalc,'go');
plot(xt,usub,'r*');
plot(xt,usubcalc,'g*');

figure(6);
subplot(2,1,1);
plot(xt,usel-ucalc,'.');
axis([min(xt) max(xt) -1 1]);
subplot(2,1,2);
plot(xt,usub-usubcalc,'.');
axis([min(xt) max(xt) -1 1]);

figure(8);
colormap(gray);
imagesc(xtid,settings.dsel,scores{1,2});
axis xy

figure(9);
colormap(gray);
imagesc(xtid(xsel),settings.dsel,scores{1,2}(:,xsel));
axis xy


%%

%
std(usel-ucalc)
%
std(usub-usubcalc)





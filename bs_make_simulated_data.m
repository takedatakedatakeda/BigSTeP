function [data, onset, cpattern, spattern, signal, noise, sd_noise]...
                                            = bs_make_simulated_data(parm)
% Make simulated data
%
% -- Input
% parm : Parameters to generate simulated data with following fields
% .T : Length of data (default = 1000)
% .N : Length of spatiotemporal patterns (default = 20)
% .K : Number of spatiotemporal patterns (default = 5)
% .Nonset : Number of onsets for each spatiotemporal pattern (default = 25)
% .minIOI : Minimum inter-onset interval of a pattern (default = 1)
% .CH : Number of channels (default = 10)
% .SNR : Signal to noise ratio (default = -5)
% .dev : Deviation of subject-specific spatiotemporal patterns from common ones (default = 0)
%
% -- Output
% data : Simulated data (T x CH)
% onset : Onsets of spatiotemporal patterns (Nonset x K)
% cpattern : Common spatiotemporal patterns (N x K x CH)
% spattern : Subject-specific spatiotemporal patterns (N x K x CH)
% signal : Signal timeseries (T x CH)
% noise : Noise (T x CH) (data = signal + noise)
% sd_noise: SD of noise
%
% 2023/08/07 Yusuke Takeda

% Default parameters
T = 1000;% Length of data
N = 20;% Length of spatiotemporal patterns
K = 5;% Number of spatiotemporal patterns
Nonset = 25;% Number of onsets for each spatiotemporal pattern
minIOI = 1;% Minimum inter-onset interval
CH = 10;% Number of channels
SNR = -5;% Signal to noise ratio
dev = 0;% Deviation of subject-specific spatiotemporal patterns from common ones

% Check input parameters
if nargin > 0
    if isfield(parm, 'T')
        T = parm.T;
        fprintf('Data length is set to %6.0f. \n', T)
    end
    if isfield(parm, 'K')
        K = parm.K;
        fprintf('Number of spatiotemporal patterns is set to %2.0f. \n', K)
    end
    if isfield(parm, 'Nonset')
        Nonset = parm.Nonset;
        fprintf('Number of onsets is set to %4.0f. \n', Nonset)
    end
    if isfield(parm, 'minIOI')
        minIOI = parm.minIOI;
        fprintf('Minimum inter-onset interval is set to %3.0f. \n', minIOI)
    end
    if isfield(parm, 'CH')
        CH = parm.CH;
        fprintf('Number of channels is set to %2.0f. \n', CH)
    end
    if isfield(parm, 'SNR')
        SNR = parm.SNR;
        fprintf('SNR is set to %2.0f. \n', SNR)
    end
    if isfield(parm, 'dev')
        dev = parm.dev;
        fprintf('Deviation is set to %1.2f. \n', dev)
    end
end

% Load common spatiotemporal patterns
pattern_file = which('cpattern.mat');
load(pattern_file, 'cpattern')
spattern = cpattern;

% Add deviation
power_cpt = sum(sum(cpattern.^2, 1), 3);
power_dv = power_cpt*dev;
mag_dv = sqrt(power_dv);
width = 5;
c = ones(width, width);
for k = 1:K
    dv = filter2(c, randn(N+2*width, CH+2*width));% Smoothed noise
    dv = dv(width+1:width+N, width+1:width+CH);
    dv = dv/sqrt(sum(dv(:).^2))*mag_dv(k);
    spattern(:,k,:) = cpattern(:,k,:)+reshape(dv,N,1,CH);
end

% Make onset timeseries
onset = bs_make_random_onset(Nonset, T, N, K, minIOI);
onset_ts = bs_make_onset_timeseries(onset, T);

% Make signal
s = zeros(T+N-1, CH, K);
for ch = 1:CH
    for k = 1:K
        s(:, ch, k) = conv(onset_ts(:, k), spattern(:, k, mod(ch-1,10)+1));
    end
end
signal = sum(s(1:T, :, :),3);

% Make noise
sd_noise = (10^(-SNR/10)/N)^0.5;
noise = randn(T, CH)*sd_noise;

% Make simulated data
data = signal+noise;

% Make spatiotemporal patterns to output
cpattern1 = zeros(N, K, CH);
spattern1 = zeros(N, K, CH);
for ch = 1:CH
    cpattern1(:, :, ch) = cpattern(:, 1:K, mod(ch-1,10)+1);
    spattern1(:, :, ch) = spattern(:, 1:K, mod(ch-1,10)+1);
end
cpattern = cpattern1;
spattern = spattern1;




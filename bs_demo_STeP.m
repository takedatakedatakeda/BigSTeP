% Demo program for STeP (Takeda et al., NeuroImage 2016, 133: 251-265)
% STeP estimates repetitive spatiotemporal patterns and their onsets from
% single-subject resting-state data.
%
% 2023/08/07 Yusuke Takeda

%% Set parameters for this simulation test

clear all
close all

T = 500;% Length of simulated data
N = 20;% Length of spatiotemporal pattern
K = 5;% Number of spatiotemporal patterns
Nonset = 10;% Number of onsets for each spatiotemporal pattern
minIOI = 20;% Minimum inter-onset interval
CH = 10;% Number of channels
SNR = 0;% Signal to noise ratio

%% Make and show simulated data

% Set parameters
sim_parm.T = T;
sim_parm.N = N;
sim_parm.K = K;
sim_parm.Nonset = Nonset;
sim_parm.minIOI =  minIOI;
sim_parm.CH = CH;
sim_parm.SNR = SNR;

% Make simulated data
[data, onset, ~, pattern, signal, noise] = bs_make_simulated_data(sim_parm);

% Show spatiotemporal patterns and their onset timeseries
figure(1);clf
ma = max(pattern(:));
mi = min(pattern(:));
onset_ts = bs_make_onset_timeseries(onset, T);
for k = 1:K
    subplot(K, 4, 4*(k-1)+1)
    imagesc(squeeze(pattern(:, k, :))')
    caxis([mi ma])
    colorbar
    if k == 1
        title('True spatiotemporal pattern')
    end
    if k == K
        xlabel('Time')
    end
    if k == fix(K/2)+1
        ylabel('Channel')
    end
    subplot(K, 4, 4*(k-1)+2:4*k)
    plot(onset_ts(:, k))
    if k == 1
        title('True onset timeseries')
    end
    if k == K
        xlabel('Time')
    end
end

% Show signal, noise, and simulated data
figure(2);clf
mi = min(data(:));
ma = max(data(:));
subplot(3, 1, 1)
imagesc(signal')
caxis([mi ma])
colorbar
title('Signal')
subplot(3, 1, 2)
imagesc(noise')
caxis([mi ma])
colorbar
title('Noise')
ylabel('Channel')
subplot(3, 1, 3)
imagesc(data')
caxis([mi ma])
colorbar
title('Simulated data (= Signal + Noise)')
xlabel('Time')


%% Estimate spatiotemporal patterns and their onsets using STeP

% Apply STeP to simulated data
STeP_parm.minIOI = minIOI;
tic
[e_onset, e_pattern] = bs_STeP(data, N, K, STeP_parm);
toc 

% Show estimation result
figure(3);clf
ma = max(e_pattern(:));
mi = min(e_pattern(:));
e_onset_ts = bs_make_onset_timeseries(e_onset, T);
for k = 1:K
    subplot(K, 4, 4*(k-1)+1)
    imagesc(squeeze(e_pattern(:, k, :))')
    caxis([mi ma])
    colorbar
    if k == 1
        title('Estimated spatiotemporal pattern')
    end
    if k == K
        xlabel('Time')
    end
    if k == fix(K/2)+1
        ylabel('Channel')
    end
    subplot(K, 4, 4*(k-1)+2:4*k)
    plot(e_onset_ts(:, k))
    if k == 1
        title('Estimated onset timeseries')
    end
    if k == K
        xlabel('Time')
    end
end

%% Adjust estimated onsets to match true ones

% Adjust estimated onsets
[a_onset, a_pattern] = bs_adjust_onset_to_ref(data, pattern, e_onset, N);

% Show adjusted results
figure(4);clf
ma = max(a_pattern(:));
mi = min(a_pattern(:));
a_onset_ts = bs_make_onset_timeseries(a_onset, T);
for k = 1:K
    subplot(K, 4, 4*(k-1)+1)
    imagesc(squeeze(a_pattern(:, k, :))')
    caxis([mi ma])
    colorbar
    if k == 1
        title('Adjusted spatiotemporal pattern')
    end
    if k == K
        xlabel('Time')
    end
    if k == fix(K/2)+1
        ylabel('Channel')
    end
    subplot(K, 4, 4*(k-1)+2:4*k)
    plot(a_onset_ts(:, k))
    if k == 1
        title('Adjusted onset timeseries')
    end
    if k == K
        xlabel('Time')
    end
end

%% Quantify estimation accuracy

% Correlation coefficient between true and estimated spatiotemporal patterns
r = bs_accuracy_of_pattern(pattern, a_pattern);
fprintf('Correlation coefficient between true and estiated patterns is %1.2f.\n', r)

% Normalized distance from true onsets
nd = bs_calc_normalized_dist(onset, a_onset, T);
fprintf('Normalized distance from true onset is %1.4f.\n', nd)

% Normalized number of estimated onsets
nn = bs_calc_normalized_num(onset, a_onset);
fprintf('Normalized number of estimated onsets is %1.2f.\n', nn)

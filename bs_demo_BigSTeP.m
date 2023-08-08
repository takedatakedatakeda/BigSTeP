% Demo program for BigSTeP (Takeda et al., NeuroImage 2019, 203: 116182)
% BigSteP estimates repetitive spatiotemporal patterns and their onsets 
% from multi-subject resting-state data.
%
% 2023/08/07 Yusuke Takeda

%% Set parameters for this simulation test

clear all
close all

Nsub = 10;% Number of subjects
T = 500;% Length of simulated data
N = 20;% Length of spatiotemporal pattern
K = 5;% Number of spatiotemporal patterns
Nonset = 10;% Number of onsets for each spatiotemporal pattern
minIOI = 20;% Minimum inter-onset interval
CH = 10;% Number of channels
SNR = 0;% Signal to noise ratio
dev = 0.1;% Deviation of subject-specific spatiotemporal patterns from common ones

%% Make simulated data

% Set parameters
sim_parm.T = T;
sim_parm.N = N;
sim_parm.K = K;
sim_parm.Nonset = Nonset;
sim_parm.minIOI =  minIOI;
sim_parm.CH = CH;
sim_parm.SNR = SNR;
sim_parm.dev = dev;

% Make simulated data
data = cell(1, Nsub);
onset = cell(1, Nsub);
spat = cell(1, Nsub);
for sub = 1:Nsub
    [data{sub}, onset{sub}, cpat, spat{sub}] = bs_make_simulated_data(sim_parm);
end

% Obtain range of values to show
ma = max(cpat(:));
mi = min(cpat(:));
for sub = 1:Nsub
    ma = max([max(spat{sub}(:)) ma]);
    mi = min([min(spat{sub}(:)) mi]);
end

% Show true spatiotemporal patterns
figure(1);clf
for k = 1:K
    subplot(K, 3, 3*(k-1)+1)
    imagesc(squeeze(cpat(:, k, :))')
    caxis([mi ma])
    colorbar
    if k == 1
        title('True common pattern')
    end
    if k == fix(K/2)+1
        ylabel('Channel')
    end
    
    for sub = 1:2
        subplot(K, 3, 3*(k-1)+1+sub)
        imagesc(squeeze(spat{sub}(:, k, :))')
        caxis([mi ma])
        colorbar
        if k == 1
            title(['True subj.' num2str(sub) '-specific pattern'])
        end
        if sub == 1 && k == K
            xlabel('Time')
        end
    end
end

% Show simulated data
figure(2);clf
for sub = 1:5
    subplot(5, 1, sub)
    imagesc(data{sub}')
    colorbar
    title(['Simulated data (Subj. ' num2str(sub) ')'])
    if sub == 5
        xlabel('Time')
    end
    if sub == 3
        ylabel('Channel')
    end
end

% Quantify how frequently patterns are ovarlapping in simulated data
p_overlap = bs_proportion_of_overlap(T, onset, N);
fprintf('Proportion of overlap is %1.2f.\n', p_overlap)

%% Apply BigSTeP to simulated data

% Estimate common and subject-specific spatiotemporal patterns
STeP_parm.minIOI = minIOI;
[e_onset, e_cpat, e_spat, e_onset_step] = bs_BigSTeP(data, N, K, STeP_parm);

% Adjust estimated onsets to true ones
[a_onset, a_cpat] = bs_adjust_onset_to_ref(data, cpat, e_onset, N);

% Estimate subject-specific spatiotemporal patterns using adjusted onsets
a_spat = cell(1, Nsub);
for sub = 1:Nsub
    a_spat{sub} = bs_est_pattern(data{sub}, a_onset{sub}, N);
end

% Obtain range of values to show
ma = max(a_cpat(:));
mi = min(a_cpat(:));
for sub = 1:Nsub
    ma = max([max(a_spat{sub}(:)) ma]);
    mi = min([min(a_spat{sub}(:)) mi]);
end

% Show estimated spatiotemporal patterns
figure(3);clf
for k = 1:K
    subplot(K, 3, 3*(k-1)+1)
    imagesc(squeeze(a_cpat(:, k, :))')
    caxis([mi ma])
    colorbar
    if k == 1
        title('Estimated common pattern')
    end
    if k == fix(K/2)+1
        ylabel('Channel')
    end
    
    for sub = 1:2
        subplot(K, 3, 3*(k-1)+1+sub)
        imagesc(squeeze(a_spat{sub}(:, k, :))')
        caxis([mi ma])
        colorbar
        if k == 1
            title(['Estimated subj.' num2str(sub) '-specific pattern'])
        end
        if sub == 1 && k == K
            xlabel('Time')
        end
    end
end

%% Extract significantly large activities in common spatiotemporal patterns

% Estimate FDR for estimated common spatiotemporal patterns
q = bs_statistical_test(data, a_onset, N);

% Extract significant activities
sig_cpat = a_cpat.*(q < 0.01);

% Show significant activities
figure(4);clf
for k = 1:K
    subplot(K, 3, 3*(k-1)+2)
    imagesc(squeeze(sig_cpat(:, k, :))')
    caxis([mi ma])
    colorbar
    if k == 1
        title('Significant activity in estimated common pattern')
    end
    if k == fix(K/2)+1
        ylabel('Channel')
    end
    if k == K
        xlabel('Time')
    end
end

%% Quantify estimation accuracy

% Correlation coefficient between true and estimated spatiotemporal patterns
r = bs_accuracy_of_pattern(spat, a_spat);
fprintf('Correlation coefficient is %1.2f.\n', r)

% Correlation coefficient of STeP result (for comparison)
a_pat_step = cell(1, Nsub);
for sub = 1:Nsub
    [~, a_pat_step{sub}] = bs_adjust_onset_to_ref(data{sub}, spat{sub}, e_onset_step{sub}, N);
end
r_step = bs_accuracy_of_pattern(spat, a_pat_step);
fprintf('Correlation coefficient (STeP) is %1.2f.\n', r_step)

% Normalized distance from true onsets
nd = bs_calc_normalized_dist(onset, a_onset, T);
fprintf('Normalized distance from true onset is %1.4f.\n', nd)

% Normalized number of estimated onsets
nn = bs_calc_normalized_num(onset, a_onset);
fprintf('Normalized number of estimated onsets is %1.2f.\n', nn)

%% Calculate contribution ratio

% Calculate contribution ratio
[contribution_ratio, label] = bs_contribution_ratio(data, a_onset, N);

% Show contribution ratio
figure(6);clf
pie(contribution_ratio, label)
title('Contribution ratio')

%% Calculate cross-correlograms from estimated onsets

% Calculate cross-correlograms
width = 10;
[cc, time_label] = bs_cross_correlogram(onset, width);

% Show cross-correlograms
figure(7);clf
for k1 = 1:K
    for k2 = 1:K
        subplot(K, K, K*(k1-1)+k2)
        plot(time_label, cc(:, k1, k2))
        axis([-width width 0 1])
        title([num2str(k2) ' triggering on ' num2str(k1)])
        if k1 == K
            xlabel('Time')
        end
    end
end

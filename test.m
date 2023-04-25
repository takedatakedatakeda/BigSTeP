clear all
close all

Nsub = 10;% Number of subjects
T = 500;% Length of simulated data (can be different across subjects)
N = 20;% Length of spatiotemporal pattern
K = 5;% Number of spatiotemporal patterns
Nonset = 10;% Number of onsets for each spatiotemporal pattern
CH = 10;% Number of channels
SNR = 5;% Signal to noise ratio
dev = 0;% Deviation of subject-specific spatiotemporal patterns from common ones

parm.T = T;
parm.N = N;
parm.K = K;
parm.Nonset = Nonset;
parm.CH = CH;
parm.SNR = SNR;
parm.dev = dev;

% Make simulated data
data = cell(1, Nsub);
onset = cell(1, Nsub);
spat = cell(1, Nsub);
sig = cell(1, Nsub);
noise = cell(1, Nsub);
for sub = 1:Nsub
    [data{sub}, onset{sub}, cpat, spat{sub}, sig{sub}, noise{sub}] = bs_make_simulated_data(parm);
end

sub = 1;

% Test for bs_est_pattern.m
e_spat = bs_est_pattern(sig{sub}, onset{sub}, N);
plot(spat{sub}(:)-e_spat(:))% Should be zero

e_cpat = bs_est_pattern(sig, onset, N);
plot(cpat(:)-e_cpat(:))% Should be zero

% Test for bs_predict_data.m
pdata = bs_predict_data(T, onset{sub}, spat{sub});
plot(sig{sub}(:)-pdata(:))% Should be zero

% Test for bs_residual_error.m
re = bs_residual_error(data{sub}, onset{sub}, N, spat{sub});
plot(noise{sub}(:)-re(:))% Should be zero

re = bs_residual_error(data, onset, N, spat);
plot(noise{sub}(:)-re{sub}(:))% Should be zero

% Test for bs_update_onset.m
[u_onset, u_cpat, error] = bs_update_onset(sig, onset, N);
for sub = 1:Nsub
    plot(onset{sub}(:)-u_onset{sub}(:))% Should be zero
    hold on
end
hold off
length(error)% Should be 1

% Test for bs_bottom_level.m
[u_onset, error] = bs_bottom_level(sig{sub}(N:end, :), N, N:T-N+1, onset{sub}, []);
plot(onset{sub}(:)-u_onset(:))% Should be zero

onset1 = onset{sub};
onset1(1, 1) = onset1(1, 1)+1;
onset1(1, 2) = onset1(1, 2)+1;
[u_onset, error] = bs_bottom_level(sig{sub}(N:end, :), N, N:T-N+1, onset1, []);
plot(onset{sub}(:)-u_onset(:))% Should be zero

% Test for proportion of overlap.m
p_overlap = bs_proportion_of_overlap(T, [1 1+N 1+2*N], N)% Should be 0
p_overlap = bs_proportion_of_overlap(T, [1 1+N/2 1+2*N/2], N)% Should be 0.5
p_overlap = bs_proportion_of_overlap(T, [1 1 1], N)% Should be 1

% Test for bs_make_onset_timeseries.m
onset_ts = bs_make_onset_timeseries(onset{sub}, T);
onset_ts = bs_make_onset_timeseries(onset, T);
onset_ts = bs_make_onset_timeseries(onset, ones(Nsub, 1)*T);
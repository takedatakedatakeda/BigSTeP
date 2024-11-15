function [onset, error] = bs_middle_level(data, N, onset_list, M, previous_onset, E, peak_time, minIOI)
% Middle level in STeP procedure
%
% -- Input
% data : Resting-state data. (T x CH)
% N : Length of spatiotemporal patterns
% onset_list : List of onsets (vector)
% M : Number of initial onsets for each spatiotemporal pattern
% previous_onset : Onsets estimated in previous loop at top level (Nonset x K)
% E : Number of repetition to estimate onsets
% peak_time : Peak time point
% minIOI : Minimum inter-onset interval
%
% -- Output
% onset : Estimated onsets of spatiotemporal patterns (Nonset x K)
% error : Power of residual error
%
% 2024/07/05 Yusuke Takeda

% Set values
T = size(data, 1);
Y = data(N:end, :);
K = size(previous_onset, 2);

% Make distribution of random variables for onsets
non_onset_list = 1:T;
non_onset_list(onset_list) = [];
if sum(previous_onset(:)) == 0
    distribution = bs_moving_average(sum(data.^2, 2), N);
else
    residual_error = bs_residual_error(data, previous_onset, N);
    distribution = bs_moving_average(sum(residual_error.^2, 2), N);
end
distribution(1:T-round(N/2)+1) = distribution(round(N/2):T);
distribution(non_onset_list) = 0;

% Search for onsets for E times with different initial onsets
error1 = zeros(E, 1);
onset1 = zeros(M, K, E);
fprintf('ite=       ')
for ite = 1:E % This loop can be parallelized by using "parfor".
    fprintf('\b\b\b\b\b\b\b%3.0f/%3.0f', ite, E)
    
    % Generate initial onsets
    initial_onset = zeros(M, K);
    for k = 1:K
        pon = previous_onset(:, k);
        pon = pon(pon>0);
        Npon = length(pon);
        distribution1 = distribution;
        for o = 1:Npon
            distribution1(max([pon(o)-minIOI+1 1]):min([pon(o)+minIOI-1 T])) = 0;
        end
        new_on = zeros(M-Npon, 1);
        for o = 1:M-Npon
            if sum(distribution1) == 0
                break
            else
                new_on(o) = bs_make_random_value(distribution1, 1);
                distribution1(max([new_on(o)-minIOI+1 1]):min([new_on(o)+minIOI-1 T])) = 0;
            end
        end
        on = sort([pon; new_on]);
        initial_onset(:, k) = on;
    end
    
    % Bottom level
    [onset1(:,:,ite), error1(ite,1)] = bs_bottom_level(Y, N, onset_list, initial_onset, peak_time, minIOI);
end

% Select the best onsets
[error, b] = min(error1);
onset = onset1(:, :, b);

% Add an onset if all onsets are removed
for k = 1:K
    if sum(onset(:, k)) == 0
        onset(1, k) = bs_make_random_value(distribution, 1);
    end
end

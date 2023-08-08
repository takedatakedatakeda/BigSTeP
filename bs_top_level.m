function [onset, pattern] = bs_top_level(data, N, K, onset_list, max_count, E, peak_time, minIOI)
% Top level in STeP procedure
%
% -- Input
% data : Resting-state data (T x CH)
% N : Length of spatiotemporal patterns
% K : Number of spatiotemporal patterns
% onset_list : Candidate list for onsets (vector)
% max_count : Maximum value of count
% E : Number of repetition to estimate onsets
% peak_time : Peak time point
% minIOI : Minimum inter-onset interval
%
% -- Output
% onset : Estimated onsets of spatiotemporal patterns (Nonset x K)
% pattern : Estimated spatiotemporal patterns (N x K x CH)
% 
% 2023/08/07 Yusuke Takeda

warning off

% Set values
max_M = round(length(onset_list)/N);% Maximum value of M
fprintf('Maximum number of onsets is %4.0f. \n', max_M)
count = 0;
m = 1;
previous_onset = zeros(1, K);

% Estimate spatiotemporal patterns and their onsets while gradually
% increasing the number of onsets
for M = 2:2:max_M
    
    fprintf('Count=%1.0f, M=%3.0f, ', count, M);
    
    % Middle level
    start_top = tic;
    [onset, error] = bs_middle_level(data, N, onset_list, M, previous_onset, E, peak_time, minIOI);
    fprintf(', Elapsed time=%0.2f sec.\n', toc(start_top))
    
    % Update results
    previous_onset = onset;
    onset_all(1:size(onset,1), :, m) = onset;
    error_all(m, 1) = error;
    if length(error_all) > 1
        if error_all(m) >= err_min
            count = count+1;
        else
            count = 0;
            err_min = error_all(m);
        end
    elseif length(error_all) == 1
        err_min = error_all;
    end
    
    % Exit if count exceed max_count
    if count > max_count
        break
    end
    
    m = m+1;
end
[~, b] = min(error_all);

% Make output
clear onset
for k = 1:K
    a = find(onset_all(:, k, b) > 0);
    onset(1:length(a), k) = onset_all(a, k, b);
end

pattern = bs_est_pattern(data, onset, N);

warning on
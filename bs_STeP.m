function [onset, pattern] = bs_STeP(data, N, K, parm)
% This propram performs SpatioTemporal Pattern estimation (STeP) 
% proposed in Takeda et al., NeuroImage 2016, 133: 251-265.
% STeP estimates repetitive spatiotemporal patterns and their onsets
% from single-subject resting-state data.
%
% -- Input
% data : Resting-state data (T x CH)
% N : Length of spatiotemporal patterns
% K : Number of spatiotemporal patterns
% parm : <Optional> Parameters with following fields
% .onset_list : Candidate list for onsets (default = whole time points)
% .max_count : Maximum value of count (default = 3)
% .E : Number of repetition to estimate onsets (default = 30)
% .peak_time : Peak time point in a pattern (default = 0; that is, arbitrarily determined)
% .minIOI : Minimum inter-onset interval of a pattern (default = 1 time point)
%
% -- Output
% onset : Estimated onsets of spatiotemporal patterns (Nonset x K)
% pattern : Estimated spatiotemporal patterns (N x K x CH)
%
% 2023/08/07 Yusuke Takeda

fprintf('----- STeP Start -----\n')

% Set values for top level
onset_list = N:size(data,1)-N+1;% Candidate list for onsets
max_count = 3;% Maximum value of count

% Set a value for middle level
E = 30;% Number of repetition to estimate onsets

% Set values for bottom level
peak_time = 0;% Peak time point
minIOI = 1;% Minimum inter-onset interval

% Check input parameters
if nargin > 3
    if isfield(parm, 'onset_list')
        onset_list = sort(intersect(onset_list, parm.onset_list));
        fprintf('Candidate list for onsets is set. \n')
    end
    if isfield(parm, 'max_count')
        max_count = parm.max_count;
        fprintf('Maximum value of count is set to %2.0f. \n', max_count)
    end
    if isfield(parm, 'E')
        E = parm.E;
        fprintf('Number of repetition to estimate onsets is set to %2.0f. \n', E)
    end
    if isfield(parm, 'peak_time')
        peak_time = parm.peak_time;
        fprintf('Peak time point is set to %2.0f. \n', peak_time)
    end
    if isfield(parm, 'minIOI')
        minIOI = parm.minIOI;
        fprintf('Minimum inter-onset interval is set to %3.0f. \n', minIOI)
    end
end

% Top level
start_step = tic;
[onset, pattern] = bs_top_level(data, N, K, onset_list, max_count, E, peak_time, minIOI);

fprintf('Total elapsed time of STeP = %0.2f min.\n', toc(start_step)/60)
fprintf('----- STeP Finish -----\n')





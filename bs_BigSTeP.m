function [onset, cpattern, spattern, onset_step] = bs_BigSTeP(data, N, K, STeP_parm)
% Estimate repetitive spatiotemporal patterns and their onsets from many subjects' resting-state data
%
% -- Input
% data : Resting-state data (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
% K : Number of spatiotemporal patterns
% STeP_parm : <Optional> STeP's parameters with following fields
% .onset_list : Candidate list for onsets
% .max_count : Maximum value of count
% .E : Number of repetition to estimate onsets
% .peak_time : Peak time point
%
% -- Output
% onset : Estimated onsets of spatiotemporal patterns (1 x Nsub cell array)
% cpattern : Estimated common spatiotemporal patterns (N x K x CH)
% spattern : Estimated subject-specific spatiotemporal patterns (1 x Nsub cell array)
% onset_step : Onsets estimated by STeP (1 x Nsub cell array)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

start_bigstep = tic;
fprintf('---------- BigSTeP Start ----------\n')

% Convert matrix to cell if data is not cell
if ~iscell(data)
    data = {data};
end

% Stage 1: Apply STeP for each subject
fprintf('Stage 1: Apply STeP for each subject\n')
if ~exist('STeP_parm', 'var')
    STeP_parm = [];
end
Nsub = length(data);
onset_step = cell(1, Nsub);
for sub = 1:Nsub
    disp(['Apply STeP to ' num2str(sub) '-th of ' num2str(Nsub) ' subjects'])
    onset_step{sub} = bs_STeP(data{sub}, N, K, STeP_parm);
end

% Stage 2: Estimate common spatiotemporal patterns
fprintf('Stage 2: Estimate common spatiotemporal patterns\n')
m_onset = bs_match_onset(data, onset_step, N);
[onset, cpattern] = bs_update_onset(data, m_onset, N);

% Stage 3: Estimate subject-specific spatiotemporal patterns
fprintf('Stage 3: Estimate subject-specific spatiotemporal patterns\n')
spattern = cell(1, Nsub);
for sub = 1:Nsub
    spattern{sub} = bs_est_pattern(data{sub}, onset{sub}, N);
end

fprintf('Total elapsed time of BigSTeP = %0.2f min.\n', toc(start_bigstep)/60)
fprintf('---------- BigSTeP Finish ----------\n')





function [q, p] = bs_statistical_test(data, onset, N, Nshuffle)
% Statistical test for estimated common spatiotemporal patterns
%
% -- Input
% data: Data (T x CH) or (1 x Nsub cell array)
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
% Nshuffle : Number of shuffling
%
% -- Output
% q : q-value (N x K x CH)
% p : p-value (N x K x CH)
%
% 2023/08/07 Yusuke Takeda

% Convert matrices to cells if not 
if ~iscell(data)
    data = {data};
    onset = {onset};
end

% Set parameters
if ~exist('Nshuffle', 'var')
    Nshuffle = 1000;
end
K = size(onset{1}, 2);
CH = size(data{1}, 2);

% Estimate patterns
pat = bs_est_pattern(data, onset, N);
apat = abs(pat);

% Shuffle IOIs
disp('Permutation test for estimated patterns')
count = pat*0;
fprintf('Permuting... %4.0f/%4.0f', 0, Nshuffle)
for shuffle = 1:Nshuffle
    fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', shuffle, Nshuffle)
    sonset = bs_shuffle_ioi(onset);
    spat = bs_est_pattern(data, sonset, N);
    aspat = abs(spat);
    count = count+(aspat >= apat);
end
fprintf('\n')

% Calculate p and q values
p = count/Nshuffle;
q = reshape(bs_convert_p2q(p(:)), N, K, CH);


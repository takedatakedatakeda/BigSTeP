function pattern = bs_est_pattern(data, onset, N)
% Estimate spatiotemporal patterns using onsets
%
% -- Input
% data : Resting-state data (T x CH) or (1 x Nsub cell array)
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
%
% -- Output
% pattern : Estimated spatiotemporal patterns (N x K x CH)
%
% 2023/08/07 Yusuke Takeda

% Make data matrix Y
[Y, onset] = bs_make_Y(data, N, onset);

% Estimate pattern matrix P
P = bs_est_P_from_Y(Y, onset, N);

% Make output
K = size(onset{1}, 2);
CH = size(Y{1}, 2);
pattern = reshape(P, N, K, CH);


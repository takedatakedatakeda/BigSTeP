function onset = bs_make_random_onset(Nonset, T, N, K)
% Make random onsets
%
% -- Input
% Nonset : Number of onsets
% T : Data length
% N : Length of spatiotemporal patterns
% K : Number of spatiotemporal patterns
%
% -- Output
% onset : Onsets of spatiotemporal patterns (Nonset x K)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

onset = zeros(Nonset, K);
onset_candidate = (N:T-N)';
Ncandidate = length(onset_candidate);

for k = 1:K
    ix = randperm(Ncandidate);
    onset(:, k) = onset_candidate(ix(1:Nonset));
end

onset = sort(onset);

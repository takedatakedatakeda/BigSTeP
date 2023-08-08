function onset = bs_make_random_onset(Nonset, T, N, K, minIOI)
% Make random onsets of spatiotemporal patterns
%
% -- Input
% Nonset : Number of onsets
% T : Data length
% N : Length of spatiotemporal patterns
% K : Number of spatiotemporal patterns
% minIOI : Minimum inter-onset interval (default = 1)
%
% -- Output
% onset : Onsets of spatiotemporal patterns (Nonset x K)
%
% 2023/08/07 Yusuke Takeda

if ~exist('minIOI', 'var')
    minIOI = 1;% Minimum inter-onset interval
end

onset = zeros(Nonset, K);
onset_candidate = (N:T-N)';

for k = 1:K
    tmp1 = onset_candidate;
    for on = 1:Nonset
        tmp2 = bs_rs(tmp1);
        onset(on, k) = tmp2(1);
        ix = find(tmp1 >= onset(on, k)-minIOI+1 & tmp1 <= onset(on, k)+minIOI-1);
        tmp1(ix) = [];
    end
end

onset = sort(onset);

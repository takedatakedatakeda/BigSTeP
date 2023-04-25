function p_overlap = bs_proportion_of_overlap(T, onset, N)
% Calculate proportion of overlap (0-1)
% (It becomes larger than zero if patterns are overlapping)
%
% - Input
%  T : Length of data (scalar) or (Nsub x 1)
%  onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
%  N : Length of spatiotemporal patterns
%
% - Output
%  p_overlap : Proportion of overlap
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Convert matrix to cell if onset is not a cell
if ~iscell(onset)
    onset = {onset};
end

% Prepare values
Nsub = length(onset);
c = ones(N, 1);
K = size(onset{1}, 2);
if length(T) == 1
    T = repmat(T, Nsub, 1);
end

% Calculate proportion of overlap
p_overlap = zeros(Nsub, 1);
for sub = 1:Nsub
    u = zeros(T(sub), 1);
    for k = 1:K
        a = find(onset{sub}(:, k) > 0);
        for on = 1:length(a)
            u(onset{sub}(a(on), k), 1) = u(onset{sub}(a(on), k), 1)+1;
        end
    end
    cover = conv(u, c);
    d_cover = sum(cover(1:T) > 0);
    d_unoverlap = N*sum(u);
    p_overlap(sub) = (d_unoverlap-d_cover)/(d_unoverlap-N);
end
p_overlap = mean(p_overlap);

function onset_ts = bs_make_onset_timeseries(onset, T)
% Make onset timeseries using onsets
%
% -- Input
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% T : Length of data (Nsub x 1)
%
% -- Output
% onset_ts : Onset timeseries (T x K) or (1 x Nsub cell array)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Convert matrix to cell if onset is not cell
c = 1;
if ~iscell(onset)
    onset = {onset};
    c = 0;
end

% Prepare values
Nsub = length(onset);
if length(T) == 1
    T = repmat(T, Nsub, 1);
end
K = size(onset{1}, 2);

% Make onset timeseries
onset_ts = cell(1, Nsub);
for sub = 1:Nsub
    onset_ts{sub} = zeros(T(sub), K);
    for k = 1:K
        a = find(onset{sub}(:, k) > 0);
        for on = 1:length(a)
            onset_ts{sub}(onset{sub}(a(on), k), k) = onset_ts{sub}(onset{sub}(a(on), k), k)+1;
        end
    end
end

% Make output
if c == 0
    onset_ts = onset_ts{1};
end
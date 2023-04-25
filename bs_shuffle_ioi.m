function sonset = bs_shuffle_ioi(onset)
% Shuffle inter-onset intervals (IOI)
%
% -- Input
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
%
% -- Output
% sonset : Shuffled onset (Nonset x K) or (1 x Nsub cell array)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Convert matrix to cell if onset is not cell
c = 1;
if ~iscell(onset)
    onset = {onset};
    c = 0;
end

% Shuffle IOI
sonset = onset;
for sub = 1:length(onset)
    sonset{sub} = in_each_sub(onset{sub});
end

% Convert cell to matrix if necessary
if c == 0
    sonset = sonset{1};
end

% Inner function to shuffle IOI
function sonset = in_each_sub(onset)

onset = sort(onset);
K = size(onset, 2);

sonset = onset*0;
for k = 1:K
    ix = find(onset(:, k) > 0);
    ioi = diff([0; onset(ix, k)]);
    sioi = in_rs(ioi);
    sonset(1:length(ix), k) = cumsum(sioi);
end
sonset = sort(sonset);

% Inner function to randomly shuffle x
function y = in_rs(x)

tmp = [rand(size(x, 1), 1), x];
tmp = sortrows(tmp, 1);
y = tmp(:, 2:end);

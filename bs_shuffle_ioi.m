function sonset = bs_shuffle_ioi(onset)
% Make surrogate onsets by shuffling inter-onset intervals (IOI)
%
% -- Input
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
%
% -- Output
% sonset : IOI-shuffled onsets (Nonset x K) or (1 x Nsub cell array)
%
% 2023/08/07 Yusuke Takeda

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
    sioi = bs_rs(ioi);
    sonset(1:length(ix), k) = cumsum(sioi);
end
sonset = sort(sonset);



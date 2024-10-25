function sonset = bs_shuffle_ioi(onset, T)
% Make surrogate onsets by shuffling inter-onset intervals (IOI)
%
% -- Input
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
%
% -- Output
% sonset : IOI-shuffled onsets (Nonset x K) or (1 x Nsub cell array)
%
% 2024/10/10 Yusuke Takeda

% Convert matrix to cell if onset is not cell
c = 1;
if ~iscell(onset)
    onset = {onset};
    c = 0;
end

% Shuffle IOI
sonset = onset;
for sub = 1:length(onset)
    sonset{sub} = in_each_sub(onset{sub}, T);
end

% Convert cell to matrix if necessary
if c == 0
    sonset = sonset{1};
end

% Inner function to shuffle IOI
function sonset = in_each_sub(onset, T)

K = size(onset, 2);

sonset = onset*0;
for k = 1:K
    ix = find(onset(:, k) > 0);
    tmp = unique([0; onset(ix, k); T]);
    ioi = diff(tmp);
    sioi = bs_rs(ioi);
    sonset(1:length(ix), k) = cumsum(sioi(1:end-1));
end
sonset = sort(sonset);



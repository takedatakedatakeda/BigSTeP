function sonset = bs_shuffle_ioi(onset, Ts)
% Make surrogate onsets by shuffling inter-onset intervals (IOI)
%
% -- Input
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% Ts: Lengths of data (scalar) of (Nsub x 1)
%
% -- Output
% sonset : IOI-shuffled onsets (Nonset x K) or (1 x Nsub cell array)
%
% 2024/11/5 Yusuke Takeda

% Convert matrix to cell if onset is not cell
c = 1;
if ~iscell(onset)
    onset = {onset};
    c = 0;
end

if length(Ts) == 1
    Ts = repmat(Ts, length(onset), 1);
end

% Shuffle IOI
sonset = onset;
for sub = 1:length(onset)
    sonset{sub} = bs_shuffle_ioi_each_sub(onset{sub}, Ts(sub));
end

% Convert cell to matrix if necessary
if c == 0
    sonset = sonset{1};
end





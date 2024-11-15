function sonset = bs_shuffle_ioi_each_sub(onset, T)
% Shuffle inter-onset intervals of a subject

K = size(onset, 2);

onset = sort(onset);
sonset = onset*0;
for k = 1:K
    ix = find(onset(:, k) > 0);
    if ix
        tmp = onset(ix, k);
        ioi = diff(tmp);
        sioi = [tmp(1); bs_rs(ioi)];
        sonset(1:length(ix), k) = cumsum(sioi);
    end
end
sonset = sort(sonset);
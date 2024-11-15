function cc = bs_calc_ccg_from_onset_ts(onset_ts, width)
% Calculate cross-correlogram from onset timeseries of spatiotemporal patterns
%
% -- Input
% onset_ts : Onsets timeseries of spatiotemporal patterns (Nt x K)
% width : width of cross-correlogram
%
% -- Output
% cc : Cross-correlogram (time x target x reference) (2*width+1 x K x K)
%
% 2024/11/06 Yusuke Takeda

[T, K] = size(onset_ts);

cc = zeros(width*2+1, K, K);
for ref = 1:K
    a = find(onset_ts(:, ref) == 1);
    a = a(a>width & a<T-width);
    Non = length(a);
    if Non
        for on = 1:Non
            cc(:, :, ref) = cc(:, :, ref)+onset_ts(a(on)-width:a(on)+width, :)/Non;
        end
    end
end

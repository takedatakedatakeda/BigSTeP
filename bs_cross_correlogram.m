function [cc, time_label] = bs_cross_correlogram(onset, width)
% Calculate cross-correlogram of onsets
%
% -- Input
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% width : width of cross-correlogram
%
% -- Output
% cc : Cross-correlogram (2*width+1 x K x K)
% time_label : Label of time (2*width+1 x 1)
%
% 2023/08/07 Yusuke Takeda

% Convert matrix to cell if onset is not cell
if ~iscell(onset)
    onset = {onset};
end

% Calculate cross-correlogram
Nsub = length(onset);
cc = 0;
for sub = 1:Nsub
    cc_each = in_each_sub(onset{sub}, width);
    cc = cc + cc_each/Nsub;
end

% Make time label
time_label = (-width:width)';

% Inner function to calculate cross-correlogram for each subject
function cc = in_each_sub(onset, width)

T = max([onset(:);0])+width;
onset_ts = bs_make_onset_timeseries(onset, T);
K = size(onset, 2);

cc = zeros(width*2+1, K, K);
for k1 = 1:K% Trigger
    for k2 = 1:K% Target
        a = find(onset(:, k1) > width);
        Nonset = length(a);
        if Nonset>0
            cc1 = zeros(width*2+1, Nonset);
            for on = 1:Nonset
                cc1(:, on) = onset_ts(onset(a(on), k1)-width:onset(a(on), k1)+width, k2);
            end
            cc(:, k1, k2) = mean(cc1, 2);
        end
    end
end


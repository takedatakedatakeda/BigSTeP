function [nn, nn_each_sub] = bs_calc_normalized_num(t_onset, e_onset)
% Calculate normalized number of onsets
%
% -- Input
% t_onset : True onset (Nonset x K) or (1 x Nsub cell array)
% e_onset : Estimated onset (Nonset x K) or (1 x Nsub cell array)
%
% -- Output
% nn : Normalized number of onsets averaged across subjects
% nn_each_sub : Normalized number of onsets for each subject
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Convert matrices to cells if onsets are not cells
if ~iscell(t_onset)
    t_onset = {t_onset};
    e_onset = {e_onset};
end

% Calculate normalized number of onsets
Nsub = length(t_onset);
nn_each_sub = zeros(Nsub, 1);
for sub = 1:Nsub
    nn_each_sub(sub) = in_each_sub(t_onset{sub}, e_onset{sub});
end
nn = mean(nn_each_sub);

% Inner function to calculate normalized number of onsets for each subject
function nn = in_each_sub(t_onset, e_onset)

K = size(t_onset, 2);

nn1 = zeros(K, 1);
for k = 1:K
    t = find(t_onset(:, k) > 0);
    a = find(e_onset(:, k) > 0);
    nn1(k) = length(a)/length(t);
end
nn = mean(nn1);
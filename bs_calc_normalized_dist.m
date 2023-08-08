function [nd, nd_each_sub] = bs_calc_normalized_dist(t_onset, e_onset, T)
% Calculate normalized distance of estimated onsets to true ones
%
% -- Input
% t_onset : True onset (Nonset x K) or (1 x Nsub cell array)
% e_onset : Estimated onset (Nonset x K) or (1 x Nsub cell array)
% T : Data length (Nsub x 1)
%
% -- Output
% nd : Normalized distance averaged across subjects
% nd_each_sub : Normalized distance for each subject
%
% 2023/08/07 Yusuke Takeda

% Convert matrices to cells if onsets are not cells
if ~iscell(t_onset)
    t_onset = {t_onset};
    e_onset = {e_onset};
end

% Prepare values
Nsub = length(t_onset);
if length(T) == 1
    T = ones(Nsub, 1)*T;
end

% Calculate normalized distance
nd_each_sub = zeros(Nsub, 1);
for sub = 1:Nsub
    nd_each_sub(sub) = in_each_sub(t_onset{sub}, e_onset{sub}, T(sub));
end
nd = mean(nd_each_sub);

% Inner function to calculate normalized distance for each subject
function nd = in_each_sub(t_onset, e_onset, T)

K = size(t_onset, 2);

nd1 = zeros(K, 1);
for k = 1:K
    e_onset1 = e_onset(e_onset(:,k)>0, k);
    ioi = T/(length(e_onset1)+1);
    t_onset1 = t_onset(t_onset(:,k)>0, k);
    Nonset = length(t_onset1);
    d1 = zeros(Nonset, 1);
    for on = 1:Nonset
        d1(on) = min(abs(t_onset1(on)-e_onset(:, k)));
    end
    nd1(k) = mean(d1)/ioi;
end
nd = mean(nd1);
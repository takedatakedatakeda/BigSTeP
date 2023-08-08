function [a_onset, a_pattern, order] = bs_adjust_onset_to_ref(data, r_pattern, onset, N, width)
% Adjust onsets to match reference spatiotemporal patterns
% (Number of reference patterns should be <= the number of patterns to be adjusted)
%
% -- Input
% data : Resting-state data (T x CH) or (1 x Nsub cell array)
% r_pattern : Reference spatiotemporal patterns (N x K x CH)
% onset : Target onsets to be adjusted (Nonset x K) or (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
% width : Width of time-shift for adjustment (default = N/2)
%
% -- Output
% a_onset : Adjusted onsets (Nonset x Npattern) or (1 x Nsub cell array)
% a_pattern : Adjusted spatiotemporal patterns (N x K x CH)
% order : Order of adjusted spatiotemporal patterns
%
% 2023/08/07 Yusuke Takeda

% Prepare data and values
if iscell(onset)
    c = 1;
else
    c = 0;
end
[Y, onset] = bs_make_Y(data, N, onset);
tK = size(onset{1}, 2);% Number of target patterns
rK = size(r_pattern, 2);% Number of reference patterns
CH = size(Y{1}, 2);
Nsub = length(Y);

if ~exist('width', 'var')
    width = fix(N/2);
end
N_list = -width:width;

rpattern1 = reshape(permute(r_pattern, [1 3 2]), N*CH, rK);% (N x CH) x rK

% Calculate distances between reference and target spatiotemporal patterns
re2 = zeros(rK, tK);
nd = zeros(rK, tK);
for k = 1:tK
    re = zeros(length(N_list), rK);
    for n = 1:length(N_list)
        onset1 = onset;
        for sub = 1:Nsub
            onset1{sub}(onset1{sub}(:,k)>0, k) = onset1{sub}(onset1{sub}(:,k)>0, k)+N_list(n);
        end
        P = bs_est_P_from_Y(Y, onset1, N);
        pattern = reshape(P, N, tK, CH);
        pattern1 = pattern(:, k, :);
        pattern1 = repmat(pattern1(:), 1, rK);
        re(n, :) = mean((rpattern1-pattern1).^2, 1);
    end
    [a, b] = min(re, [], 1);
    re2(:, k) = a';
    nd(:, k) = b';
end

% Obtain order of spatiotemoporal patterns
mK = min([rK tK]);
jun = perms(1:tK);
s = zeros(size(jun, 1), 1);
for j = 1:size(jun, 1)
    for k = 1:mK
        s(j) = s(j)+re2(k, jun(j,k));
    end
end
[~, b] = min(s);
order = jun(b, :);

% Shift onsets along time
a_onset = onset;
for sub = 1:Nsub
    a_onset{sub} = onset{sub}(:,order);
    for k = 1:mK
        a_onset{sub}(a_onset{sub}(:, k)>0, k) = ...
            a_onset{sub}(a_onset{sub}(:, k)>0, k)+N_list(nd(k, order(k)));
    end
end

% Obtain adjusted spatiotemporal patterns
if nargout > 1
    P = bs_est_P_from_Y(Y, a_onset, N);
    a_pattern = reshape(P, N, tK, CH);
end

if c == 0
    a_onset = a_onset{1};
end




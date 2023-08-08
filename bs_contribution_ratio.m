function [contribution_ratio, label] = bs_contribution_ratio(data, onset, N)
% Calculate contribution ratio
%
% -- Input
% data : Resting-state data (T x CH) or (1 x Nsub cell array)
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
%
% -- Output
% contribution_ratio : Contribution Ratio (1 x K+1)
% label : Label (1 x K+1)
%
% 2023/08/07 Yusuke Takeda

% Convert matrices to cells if data is not cell
if ~iscell(data)
    data = {data};
    onset = {onset};
end

% Set parameters
Nsub = length(data);
K = size(onset{1}, 2);

% Estimate subject-specific spatiotemporal patterns
pattern = cell(1, Nsub);
for sub = 1:Nsub
    pattern{sub} = bs_est_pattern(data{sub}, onset{sub}, N);
end

% Calculate contribution ratio of all the patterns
power_orig = 0;
for sub = 1:Nsub
    power_orig = power_orig+sum(data{sub}(:).^2);
end
residual_error = bs_residual_error(data, onset, N, pattern);
power_resi = 0;
for sub = 1:Nsub
    power_resi = power_resi+sum(residual_error{sub}(:).^2);
end

contribution_ratio_all = (power_orig-power_resi)/power_orig;

% Calculate contribution ratio for each pattern
power_resi_pat = zeros(1, K);
for k = 1:K
    onset1 = cell(1, Nsub);
    pattern1 = cell(1, Nsub);
    for sub = 1:Nsub
        onset1{sub} = onset{sub}(:, k);
        pattern1{sub} = pattern{sub}(:, k, :);
    end
    re = bs_residual_error(data, onset1, N, pattern1);
    for sub = 1:Nsub
        power_resi_pat(k) = power_resi_pat(k)+sum(re{sub}(:).^2);
    end
end

contribution = power_orig-power_resi_pat;
contribution_ratio = contribution/sum(contribution)*contribution_ratio_all;
contribution_ratio = [contribution_ratio, 1-sum(contribution_ratio)];

% Make label
label = cell(1, K+1);
for k = 1:K
    label{k} = ['Pattern ' num2str(k)];
end
label{k+1} = 'Residual';



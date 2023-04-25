function [m_onset, m_pattern] = bs_match_onset(data, onset, N, width, Nref)
% Match onsets across subjects
%
% -- Input
% data : Multi-subject resting-state data (1 x Nsub cell array)
% onset : Onsets of spatiotemporal patterns (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
% width : Width of time-shift for adjustment
% Nref: Number of subjects from which reference subject is selected
%
% -- Output
% m_onset : Matched onsets (1 x Nsub cell array)
% m_pattern : Matched spatiotemporal patterns (1 x Nsub cell array)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

disp('Adjust onsets across subjects')
warning off

% Set values
Nsub = length(data);
if ~exist('witdh', 'var') || isempty(width)
    width = fix(N/2);
end
if ~exist('Nref', 'var')
    Nref = Nsub;
end

% Select reference subject from Nref subjects
ref_list = round(linspace(1, Nsub, Nref));
d = zeros(Nref, Nref);
for ref1 = 1:Nref
    r_pattern = bs_est_pattern(data{ref_list(ref1)}, onset{ref_list(ref1)}, N);
    for ref2 = 1:Nref
        [~, a_pattern] = bs_adjust_onset_to_ref(data{ref_list(ref2)}, r_pattern, onset{ref_list(ref2)}, N, width);
        d(ref1, ref2) = sum((r_pattern(:)-a_pattern(:)).^2);
    end
end
[~, b] = min(mean(d, 2));
rsub = ref_list(b);

% Calculate reference patterns
r_pattern = bs_est_pattern(data{rsub}, onset{rsub}, N);

% Adjust onsets to match reference patterns
m_onset = cell(1, Nsub);
m_pattern = cell(1, Nsub);
for sub = 1:Nsub
    [m_onset{sub}, m_pattern{sub}] = bs_adjust_onset_to_ref(data{sub}, r_pattern, onset{sub}, N, width);
end

warning on
function [r, r_each_sub] = bs_accuracy_of_pattern(t_pattern, e_pattern)
% Evaluate accuracy of estimated spatiotemporal patterns by correlation coefficient
%
% -- Input
% t_pattern : True spatiotemporal patterns (N x K x CH) or (1 x Nsub cell array)
% e_pattern : Estimated spatiotemporal patterns (N x K x CH) or (1 x Nsub cell array)
%
% -- Output
% r : Correlation coefficient averaged across subjects
% r_each_sub : Correlation coefficient for each subject (Nsub x 1)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Convert matrices to cells if they are not cells
if ~iscell(t_pattern)
    t_pattern = {t_pattern};
    e_pattern = {e_pattern};
end

% Calculate correlation coefficient
Nsub = length(t_pattern);
r_each_sub = zeros(Nsub, 1);
for sub = 1:Nsub
    r_each_sub(sub) = bs_cc(t_pattern{sub}(:), e_pattern{sub}(:));
end

% Average correlation coefficients across subjects
r = mean(r_each_sub);

function [P, U] = bs_est_P_from_Y(Y, onset, N)
% 2023/08/07 Yusuke Takeda

% Convert matrices to cells if Y is not cell
if ~iscell(Y)
    Y = {Y};
    onset = {onset};
end

% Set values
Nsub = length(Y);
K = size(onset{1}, 2);
CH = size(Y{1}, 2);

% Make U'U and U'Y
U = cell(1, Nsub);
UU = zeros(N*K, N*K);
UY = zeros(N*K, CH);
for sub = 1:Nsub
    U{sub} = bs_make_U(Y{sub}, onset{sub}, N);
    UU = UU+U{sub}'*U{sub};
    UY = UY+U{sub}'*double(Y{sub});
end

% Solve Y=UA
P = UU\UY;
P(isnan(P)) = 0;
P(isinf(P)) = 0;

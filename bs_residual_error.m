function residual_error = bs_residual_error(data, onset, N, pattern)
% Calculate residual error
%
% -- Input
% data : Resting-state data (T x CH) or (1 x Nsub cell array)
% onset : Onsets of spatiotemporal patterns (Nonset x K) or (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
% pattern : Subject-specific spatiotemporal patterns (T x CH) or (1 x Nsub cell array)
%
% -- Output
% residual_error : Residual error (T x CH) or (1 x Nsub cell array)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Convert matrices to cells if they are not cells
c = 1;
if ~iscell(data)
    data = {data};
    onset = {onset};
    if exist('pattern', 'var')
        pattern = {pattern};
    end
    c = 0;
end

% Make data matrix Y
[Y, onset] = bs_make_Y(data, N, onset);

% Estimate pattern matrix P
Nsub = length(Y);
P = cell(1, Nsub);
U = cell(1, Nsub);
if exist('pattern', 'var')
    [~, K, CH] = size(pattern{1});
    for sub = 1:Nsub
        U{sub} = bs_make_U(Y{sub}, onset{sub}, N);
        P{sub} = reshape(pattern{sub}, N*K, CH);
    end
else
    for sub = 1:Nsub
        [P{sub}, U(sub)] = bs_est_P_from_Y(Y{sub}, onset{sub}, N);
    end
end

% Calculate residual error
residual_error = data;
for sub = 1:Nsub
    residual_error{sub}(N:end, :) = Y{sub}-U{sub}*P{sub};
end

% Convert cell to matrix if data was not a cell
if c == 0
    residual_error = residual_error{1};
end






function predicted_data = bs_predict_data(T, onset, pattern)
% Predict data from spatiotemporal patterns and their onsets
%
% -- Input
% data_length : Length of data
% onset : Onsets of spatiotemporal patterns (Nonset x K)
% pattern : Spatiotemporal patterns (N x K x CH)
%
% -- Output
% predicted_data : Predicted data (T x CH)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Set parameters
[N, K, CH] = size(pattern);

% Make predicted data
dat = zeros(T+N-1, K, CH);
for k = 1:K
    u = zeros(T, 1);
    a = find(onset(:, k) > 0);
    for tr = 1:length(a)
        u(onset(a(tr), k), 1) = u(onset(a(tr), k), 1)+1;
    end
    for ch = 1:CH
        dat(:, k, ch) = conv(u, pattern(:, k, ch));
    end
end
predicted_data = permute(sum(dat(1:T, :, :), 2), [1, 3, 2]);







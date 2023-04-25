function y = bs_moving_average(data, width)
% Calculate moving average
%
% -- Input
% data : Data (T x CH)
% width : Width of moving average
%
% -- Output
% y : Result (T x CH)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Make width odd
if mod(width, 2) == 0
    width = width+1;
end

% Calculate moving average
[T, CH] = size(data);
c = ones(width, 1)/width;
s = zeros(T+width-1, CH);
for ch = 1:CH
    s(:, ch) = conv(data(:, ch), c);
end

% Make output
y = zeros(T, CH);
y((width-1)/2+1:T-(width-1)/2, :) = s(width:T, :);


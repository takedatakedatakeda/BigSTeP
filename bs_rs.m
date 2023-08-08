function y = bs_rs(x)
% Randomly shuffle matrix x in column direction
%
% -- Input
% x : Matrix
%
% -- Output
% y : Shuffled x
%
% 2023/08/07 Yusuke Takeda

tmp = [rand(size(x, 1), 1), x];
tmp = sortrows(tmp, 1);
y = tmp(:, 2:end);
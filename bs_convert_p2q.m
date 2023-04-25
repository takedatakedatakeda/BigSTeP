function q = bs_convert_p2q(p)
% Calculate false discovery rate (FDR) by the method of Storey and Tibshrani (2003)
%
% -- Input
% p : p-value (Ntest x 1)
%
% -- Output
% q : q-value (Ntest x 1)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

m = length(p);
[sp, ix] = sort(p);

% Estimate proportion of null p-values pi0
l_list = (0:0.01:0.95)';
pi01 = l_list;
l1 = 1;
for l = l_list'
    pi01(l1) = length(sp(sp > l))/m/(1-l);
    l1 = l1+1;
end
a = polyfit(l_list, pi01, 2);
pi0 = min([1 polyval(a, 1)]);

% If pi0 < 0, estimate pi0 again by a simple method
if pi0 < 0
    l = 0.5;
    pi0 = length(sp(sp > l))/m/(1-l);
end

% Convert p-values to q-values based on pi0
q1 = pi0*sp;
for i = m-1:-1:1
    q1(i) = min([pi0*m*sp(i)/i, q1(i+1)]);
end
q = sortrows([q1, ix], 2);
q = q(:, 1);


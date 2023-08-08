function r = bs_make_random_value(distribution, Nvalue)
% Sample random values from arbitrary distribution
%
% -- Input
% distribution : Probabilistic distribution (vector)
% Nvalue : Number of random values to sample
%
% -- Output
% r : randome values (Nvalue x 1)
%
% 2023/08/07 Yusuke Takeda

b = cumsum(distribution)/sum(distribution);

r = zeros(Nvalue, 1);
for n = 1:Nvalue
    while 1
        aa = find(b > rand);
        if ~isempty(aa)
            r(n, 1) = aa(1);
            break
        end
    end
end
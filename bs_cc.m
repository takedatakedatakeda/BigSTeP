function [r, p] = bs_cc(x, y)
% Calculate correlation coefficient for each channel
%
% -- Input
% x, y : Data (Nsample x Nchannel)
%
% -- Output
% r : Correlation coefficient (1 x Nchannel)
% p : P-value (1 x Nchannel)
%
% Copyright (C) 2019, Yusuke Takeda, ATR, takeda@atr.jp

% Calculate correlation coefficients
[Nsample, CH] = size(x);
zx = (x-repmat(mean(x, 1), Nsample, 1))./repmat(std(x), Nsample, 1);
zy = (y-repmat(mean(y, 1), Nsample, 1))./repmat(std(y), Nsample, 1);
r = sum(zx.*zy, 1)./(Nsample-1);
r(isnan(r)) = 0;
r(isinf(r)) = 0;

% Calculate p-values if necessary
if nargout == 2
    fun = @(x, n) 1./(n.^0.5.*beta(n./2,1/2)).*(1+x.^2./n).^(-(n+1)./2);
    t = r.*sqrt((Nsample-2)./(1-r.^2));
    t(t<0) = -t(t<0);
    p = zeros(1, CH);
    for ch = 1:CH
        p(ch) = integral(@(x) fun(x, Nsample-2), t(ch), Inf)*2;
    end
end
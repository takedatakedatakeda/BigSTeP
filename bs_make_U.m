function U = bs_make_U(Y, onset, N)
% 2023/08/07 Yusuke Takeda

% Set values
K = size(onset, 2);
T1 = size(Y, 1);

% Make coefficient matrix U
tr1 = 0;
i1 = zeros(length(find(onset(:)))*N, 1);
j1 = i1;
s1 = i1;
for k = 1:K
    on = onset(onset(:,k)>=N & onset(:,k)<=T1, k);
    Non = length(on);
    for n = 1:N
        i1(tr1+1:tr1+Non) = on-N+n;
        j1(tr1+1:tr1+Non) = N*(k-1)+n;
        s1(tr1+1:tr1+Non) = 1;
        tr1 = tr1+Non;
    end
end
U = sparse(i1(1:tr1), j1(1:tr1), s1(1:tr1), T1, N*K);
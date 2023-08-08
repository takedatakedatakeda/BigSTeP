function [Y, onset] = bs_make_Y(data, N, onset)
% 2023/08/07 Yusuke Takeda

% Make data matrix Y
if iscell(data)
    Nsub = length(data);
    Y = cell(1, Nsub);
    for sub = 1:Nsub
        Y{sub} = data{sub}(N:end, :);
    end
else
    Y = {data(N:end, :)};
end

% Convert matrix to cell if onset is not cell
if exist('onset', 'var')
    if ~iscell(onset)
        onset = {onset};
    end
end



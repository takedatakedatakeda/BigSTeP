function [onset, cpattern, error] = bs_update_onset(data, onset, N)
% Update onsets to estimate common spatiotemporal patterns
%
% -- Input
% data : Multi-subject resting-state data (1 x Nsub cell array)
% onset : Onsets of spatiotemporal patterns (1 x Nsub cell array)
% N : Length of spatiotemporal patterns
%
% -- Output
% onset : Updated onsets (1 x Nsub cell array)
% cpattern : Common spatiotemporal patterns (N x K x CH)
% error : Amount of error
%
% 2023/08/07 Yusuke Takeda

disp('Update onsets')

% Set values
Nsub = length(data);
K = size(onset{1}, 2);
CH = size(data{1}, 2);

% Make data matrix Y
Y = bs_make_Y(data, N);

% Make onset list
onset_list = cellfun(@(x) N:size(x,1)-N+1, data, 'uniformoutput', false);

Ms = zeros(1, Nsub);
Ts = zeros(1,Nsub);
for sub = 1:Nsub
    Ts(sub) = size(data{sub}, 1);
    onset1 = onset{sub};
    for k = 1:K
        on = intersect(onset_list{sub}, onset1(:,k));
        onset1(:, k) = 0;
        onset1(1:length(on), k) = on;
    end
    onset1 = sort(onset1);
    ix = find(sum(onset1, 2) > 0);
    onset{sub} = onset1(ix, :);
    Ms(sub) = length(ix);
end
M = max(Ms);

% Alternately iterate updates of spatiotemporal patterns and a onset
ite = 1;
error(ite) = in_residual_error(Y, onset, N);
onset2 = onset;
while 1
    onset = cellfun(@sort, onset, 'uniformoutput', false);
    for k = 1:K
        ponsets = cellfun(@(x) x(1)-1, onset_list);
        for on = 1:M
            % Find sets having more onsets than on-1
            active_sub = find(Ms >= on);
            
            % Update spatiotemporal patterns and residual error using
            % onsets except for target onsets
            onset1 = onset;
            for sub = active_sub
                onset1{sub}(on, k) = 0;
            end
            [P, U] = bs_est_P_from_Y(Y, onset1, N);
            pattern = reshape(P, N, K, CH);
            pattern1 = squeeze(pattern(:, k, :));
            
            % Search for time point when residual error
            % contains target spatiotemporal pattern
            for sub = active_sub
                if onset{sub}(on, k) > 0
                    
                    R = bsxfun(@minus, Y{sub}, U{sub}*P);
                    if on == Ms(sub)
                        tlist = [0 onset_list{sub}(onset_list{sub} > ponsets(sub))];
                    else
                        tlist = [0 onset_list{sub}(onset_list{sub} > ponsets(sub) & onset_list{sub} < onset{sub}(on+1,k))];
                    end
                    r1 = R(tlist(2)-N+1:tlist(end), :);
                    pp = sum(pattern1(:).^2);
                    xy = sum(real(ifft(fft([r1;zeros(N-1,CH)]).*fft([flipud(pattern1);zeros(tlist(end)-tlist(2)+N-1,CH)]))),2);
                    e = [0; -2*xy(tlist(2:end)-tlist(2)+N)+pp];
                    [~, b] = min(e);
                    
                    % Update target onset
                    onset{sub}(on, k) = tlist(b);
                    if onset{sub}(on, k) > 0
                        ponsets(sub) = onset{sub}(on, k);
                    end
                    
                end
            end
        end
    end
    
    ite = ite+1;
    error(ite) = in_residual_error(Y, onset, N);
    if ~(error(ite) < error(ite-1))
        error = error(1:ite-1);
        onset = onset2;
        break
    else
        onset2 = onset;
    end
end

% Estimate patterns
P = bs_est_P_from_Y(Y, onset, N);
cpattern = reshape(P, N, K, CH);


% Inner function to calculate residual error of all the data
function error = in_residual_error(Y, onset, N)

% Set value
Nsub = length(Y);

% Calculate patterns and onset timeseries
[P, U] = bs_est_P_from_Y(Y, onset, N);

% Calculate residual error
error = 0;
for sub = 1:Nsub
    R = bsxfun(@minus, Y{sub}, U{sub}*P);
    error = error+sum(R(:).^2);
end

return



function [onset, error] = bs_bottom_level(Y, N, onset_list, initial_onset, peak_time, minIOI)
% Bottom level in STeP procedure
%
% -- Input
% Y : Data matrix (T-N+1 x CH)
% N : Length of spatiotemporal patterns
% onset_list : List of onsets (vector)
% initial_onset : Initial values for onsets (Nonset x K)
% peak_time : Peak time point
% minIOI : Minimum inter-onset interval
%
% -- Output
% onset : Estimated onsets (Nonset x K)
% error : Power of residual error
%
% 2023/08/07 Yusuke Takeda

% Set values
[M, K] = size(initial_onset);
onset = initial_onset;
CH = size(Y, 2);

% Alternately iterate updates of spatiotemporal patterns and a onset
ite = 1;
error(ite) = in_residual_error(Y, onset, N);
onset2 = onset;
while 1
    onset = sort(onset);
    for k = 1:K
        ponset = onset_list(1)-1;
        ix_nonzero = find(onset(:,k) > 0)';
        for on = ix_nonzero
            % Update spatiotemporal patterns and residual error using
            % onsets except for target onset
            onset1 = onset;
            onset1(on, k) = 0;
            [~, R, P] = in_residual_error(Y, onset1, N);
            pattern = reshape(P, N, K, CH);
            pattern1 = squeeze(pattern(:, k, :));
            
            % Search for time point when residual error
            % contains target spatiotemporal pattern
            if on == M
                ix = find(onset_list > ponset+minIOI-1);
            else
                ix = find(onset_list>ponset+minIOI-1 & onset_list<onset(on+1, k));
            end
            if ix
                tlist = [0 onset_list(ix)];
                r1 = R(tlist(2)-N+1:tlist(end), :);
                pp = sum(pattern1(:).^2);
                xy = sum(real(ifft(fft([r1;zeros(N-1,CH)]).*fft([flipud(pattern1); zeros(tlist(end)-tlist(2)+N-1, CH)]))), 2);
                e = [0; -2*xy(tlist(2:end)-tlist(2)+N)+pp];
                [~, b] = min(e);
                onset(on, k) = tlist(b);
            else
                onset(on, k) = 0;
            end
            if onset(on, k) > 0
                ponset = onset(on, k);
            end
        end
    end
    
    ite = ite+1;
    error(ite) = in_residual_error(Y, onset, N);
    if ~(error(ite) < error(ite-1))
        error = error(ite-1);
        onset = onset2;
        break
    else
        onset2 = onset;
    end
end

% Shift onsets
if peak_time
    [~, ~, P] = in_residual_error(Y, onset, N);
    pattern = reshape(P, N, K, CH);
    power = mean(pattern.^2, 3);
    [~, b] = max(power, [], 1);
    d = peak_time-b;
    for k = 1:K
        a = find(onset(:, k) > 0);
        onset(a, k) = onset(a, k)-d(k);
        onset(a, k) = onset(a, k).*ismember(onset(a, k),onset_list);
    end
    onset = sort(onset);
    error = in_residual_error(Y, onset, N);
end

% Inner function to calculate residual error
function [error, R, P] = in_residual_error(Y, onset, N)

% Calculate patterns and onset timeseries
[P, U] = bs_est_P_from_Y(Y, onset, N);

% Calculate residual error
R = bsxfun(@minus, Y, U{1}*P);
error = sum(R(:).^2);







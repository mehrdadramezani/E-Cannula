function [mag_true,dt_true] = obtain_mag_delay_matrix_v2(env_true,pos,mag_type)
interval_mask = 21:120; % 1-150 is 150 ms interval. 21:120 is 100 ms interval
maxlag = 15; % max lag is 15 ms
nECoG_ch = size(env_true,2);

% use the true envelop to compute the magnitude matrix
n_nonstat = size(env_true,1);
mag_true = nan(n_nonstat,16);
temp_true = mean(env_true.^2,3);  % mean power
if strcmp(mag_type, 'raw')
    % do nothing
elseif strcmp(mag_type, 'ch_minmax')
    temp_true = temp_true ./ max(temp_true,[],2); % normalize to the max channel
elseif strcmp(mag_type, 'event_mean')
    temp_true = temp_true./ mean(temp_true); % normalize to the mean over trials
else
    error('mag_type needs to be specified (raw, ch_minmax, event_minmax)'); 
end
mag_true(:,pos) = temp_true;
mag_true = reshape(mag_true, n_nonstat,4,4);

% use the true envelop to compute the delay matrix
dt_true = nan(n_nonstat,16);
for n = 1:n_nonstat
    temp_mean = mean(squeeze(env_true(n,:,:)))';
    temp1 = squeeze(env_true(n,:,:));
    [rvals, lags] = myxcorr(temp_mean, temp1, interval_mask, maxlag);
    for c = 1:nECoG_ch
        [rmax,ind] = max(rvals(c,:));
        if rmax > 0.5 % the waveform is very similar and this channel has reliable time delay estimate (corr > 0.5 uV after smoothing)
            dt_true(n,pos(c)) = lags(ind); % ms         
        else
            dt_true(n,pos(c)) = nan; % ms
        end
    end
end
dt_true = reshape(dt_true, n_nonstat,4,4); 
end
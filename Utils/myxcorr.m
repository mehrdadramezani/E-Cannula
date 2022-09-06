function [rvals,lags] = myxcorr(seq1, seq2, interval, maxlag)
% seq1 is the stationary series
% seq2 is the sliding series
% interval is for seq1 (onset:offset)
% rvals is the correlation value
% lags is the time steps that seq2 falls behind seq1

lags = (-maxlag:maxlag);
rvals = zeros(size(seq2,1),length(lags));
for r = 1:length(lags)
    offset = lags(r);
    rvals(:,r) = corr(seq1(interval),seq2(:,offset+interval)');
end
end
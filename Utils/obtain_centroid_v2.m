function [template_delays,template_mags,templates_true]=obtain_centroid_v2(clust_ind, feat_all_all,T, env_true_all)
templates_true = nan(length(clust_ind),16,T);
template_delays = nan(length(clust_ind),16);
template_mags = nan(length(clust_ind),16);
cnum = length(clust_ind);
for c = 1:cnum
    inds = clust_ind{c};
    if length(inds) > 1
        temp_true = env_true_all(inds,:,:);              % using true envelop
        templates_true(c,:,:) = mean(temp_true);
        temp_true = feat_all_all(inds,17:end); % the power features
        template_mags(c,:) = mean(temp_true);
    else
        temp_true = squeeze(env_true_all(inds,:,:));              % using true envelop
        templates_true(c,:,:) = temp_true;
        temp_true = feat_all_all(inds,17:end); % the power features
        template_mags(c,:) = temp_true;
    end

    lag_true = feat_all_all(inds,1:16); % the lag features
    template_delays(c,:) = mean(lag_true);
end
template_delays = reshape(template_delays, cnum, 4, 4);
template_mags = reshape(template_mags, cnum, 4, 4);
end
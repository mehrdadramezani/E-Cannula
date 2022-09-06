function [template_delays,template_mags,templates_true]=obtain_centroid_v2(clust_ind, feat_all_all,T, env_true_all)
% templates_recon = nan(length(clust_ind),16,T);
templates_true = nan(length(clust_ind),16,T);
template_delays = nan(length(clust_ind),16);
template_mags = nan(length(clust_ind),16);
% nstep = size(env_true_all,3);
cnum = length(clust_ind);
% nECoG_ch = size(env_true_all,2);
for c = 1:cnum
    inds = clust_ind{c};
    if length(inds) > 1
        temp_true = env_true_all(inds,:,:);              % using true envelop
%         temp_true = temp_true ./ max(reshape(temp_true,size(temp_true,1),[]),[],2);
%         temp_recon = env_recon_all(inds,:,:);            % using reconstruction
%         temp_recon = temp_recon ./ max(reshape(temp_recon,size(temp_recon,1),[]),[],2);
%         templates_recon(c,:,:) = mean(temp_recon);
        templates_true(c,:,:) = mean(temp_true);
        
%         temp_recon = mean(env_recon_all(inds,:,:).^2,3);            % using reconstruction
%         temp_recon = temp_recon ./ max(temp_recon,[],2);
        temp_true = feat_all_all(inds,17:end); % the power features
        template_mags(c,:) = mean(temp_true);
    else
        temp_true = squeeze(env_true_all(inds,:,:));              % using true envelop
%         temp_true = temp_true ./ max(temp_true(:),[],2);
%         temp_recon = squeeze(env_recon_all(inds,:,:)); % using reconstruction
%         temp_recon = temp_recon ./ max(temp_recon(:),[],2);
%         templates_recon(c,:,:) = temp_recon;
        templates_true(c,:,:) = temp_true;
        
%         temp_recon = mean(squeeze(env_recon_all(inds,:,:)).^2,2);            % using reconstruction
%         temp_recon = temp_recon ./ max(temp_recon,[],2);
        temp_true = feat_all_all(inds,17:end); % the power features
        template_mags(c,:) = temp_true;
    end
%     temp_mean = interp1(1:nstep,mean(squeeze(templates_recon(c,:,:))),linspace(1,nstep,nstep*2)); % upsample to 1 kHz
%     for n = 1:nECoG_ch
%         temp1 = squeeze(templates_true(c,n,:));
%         [r,lags] = xcorr(normalize(temp1),normalize(temp_mean),30,'unbiased'); % max lag 30 ms
%         [~,ind] = max(r);
    lag_true = feat_all_all(inds,1:16); % the lag features
    template_delays(c,:) = mean(lag_true);
%     end
end
template_delays = reshape(template_delays, cnum, 4, 4);
template_mags = reshape(template_mags, cnum, 4, 4);
end
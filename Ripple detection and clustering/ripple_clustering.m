% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for "E-Cannula reveals anatomical diversity
% in sharp-wave ripples as a driver for the recruitment of distinct
% hippocampal assemblies" published in Cell Reports.
% (C) Mehrdad Ramezani, Kuzum Lab, University of California San Diego
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code clusters the detected ripples into four distinct groups with
% different spatiotemporal profiles (Local, Global, Travelling, Stationary)
% on example LFP data.
% The result (???????) is saved to
% mat file.


%% Load the LFP data and detected SWRs
load('ripple_env.mat');
load('refined_SWRs.mat');

%% obtain the ripple onset, offset, and center time
ripple_on = ripple_result_combine.ripple_start(true_ripples)/1e3;
ripple_off = ripple_result_combine.ripple_end(true_ripples)/1e3;
ripple_centers = (ripple_on + ripple_off)/2;

%% Generate the envelop data around during each ripple
n_ripple_samples = length(ripple_centers);
ripple_centers = round(ripple_centers*fs);
ripple_on = round(ripple_on*fs);
ripple_off = round(ripple_off*fs);
% smooth the envelop data with 50 ms time window
Y_true = ripple_env'; clearvars ripple_env;
Y_true = smoothdata(Y_true, 2, 'movmean', 50);

% generate the input data
nECoG_ch = size(Y_true,1);
env_true = zeros(n_ripple_samples,nECoG_ch,0.15*fs); % 150 ms
% specify which interval to consider in later feature extraction phase.
% not enabled for current version.
mask_rp = true(n_ripple_samples,nECoG_ch,0.15*fs);
for s = 1:(n_ripple_samples)
    interval = ripple_centers(s)+(1:0.15*fs)-ceil(0.15*fs/2);
    temp = Y_true(:,interval);
    env_true(s,:,:) = temp;
    interval = (ripple_on(s):ripple_off(s)) - ripple_centers(s)+75; % mark the onset to offset intervals
    interval(interval <=0 | interval>150) = [];
    mask_rp(s,:,interval) = true;
end

%% Extract features (time delay + mean power)
% Obtain the magnitude and delay distribution matrices for non-stationary ripples
n_ripples = size(env_true,1);
mag_mode = 'ch_minmax'; % 'raw', 'ch_minmax', 'event_mean'
[mag_true,dt_true] = obtain_mag_delay_matrix_v2(env_true,pos,mag_mode);
ch_min = 12;
exclude_ind = [];
for c = 1:n_ripples
    data = squeeze(dt_true(c,:,:));
    if sum(~isnan(data(:)))<=ch_min % if the number of reliable channels is <= 12, exclude those channels
        exclude_ind = [exclude_ind,c];
    else                           % interpolate the missing channels using other data
        data = inpaint_nans(data,5);
        dt_true(c,:,:) = data;
    end
end
ind_keep = setdiff(1:n_ripples,exclude_ind);
n_ripples_keep = length(ind_keep);
env_true_keep = env_true(ind_keep,:,:);
power_feat = reshape(mag_true,n_ripples,[]);
delay_feat = reshape(dt_true,n_ripples,[]);
feat_all = [delay_feat, power_feat];

%% Separate the complex vs. simple ripples (DTW greater than threshold)
probe = normalize(exp(-((1:150)-76).^2/(2*20^2)),'range');
active_level = 3.5*std(sum(env_true_keep,3)); % the active level of each channel
n_ch = size(env_true_keep,2);
outlier_ind = [];
dist = zeros(n_ripples_keep,n_ch);
for n = 1:n_ripples_keep
    ind = squeeze(mask_rp(n,1,:));
    data = squeeze(env_true_keep(n,:,ind));
%         data = smoothdata(data','movmean',5)'; % smooth envelop at 10 ms
    for i = 1:n_ch
        dist(n,i) = dtw(max(data(i,:))*probe,data(i,:));
    end
    if any(dist(n,:)>active_level)
        outlier_ind = [outlier_ind, n];
    end
end
good_ind = setdiff(ind_keep, ind_keep(outlier_ind));
feat_good = feat_all(good_ind,:);
    
%% Perform clustering for simple ripples
feat_good_delay = feat_good(:,1:16); % the delay features
feat_good_power = feat_good(:,17:end); % the power features
feat_good_balanced = [feat_good_delay./std(feat_good_delay(:))*std(feat_good_power(:)), feat_good_power];
    

% Align the delay data to its median
feature_rescale = feat_good_balanced;
feature_rescale(:,1:16) = feat_good_balanced(:,1:16) - median(feat_good_balanced(:,1:16),2);
feat_all(:,1:16) = feat_all(:,1:16) - median(feat_all(:,1:16),2);

%% Now start the clustering
% now do the k-means clustering with different number of clusters to find optimal number of clusters
nclust_max = 30;
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);
error_val = zeros(1,nclust_max);
tic;
for k = 1:nclust_max
    [idx,C,sumd] = kmeans(feature_rescale,k,'Options',options,'Distance','cityblock','MaxIter',1000,...
        'Display','off','Replicates',6);
    error_val(k) = sum(sumd);
    toc;
end
figure('pos',[10,10,300,400]);
plot(1:nclust_max,error_val,'LineWidth',1.5);
xlabel('Number of clusters')
ylabel('Within-cluster sum of distance');
set(gca,'FontSize',12,'FontWeight','Bold');

% set the number of clusters and final k-means clustering with a specified number of clusters
rng(666);
cnum = 10; % Number of clusters (it should be optimized based on your data)
[idx_best,C,sumd] = kmeans(feature_rescale,cnum,'Options',options,'Distance','cityblock','MaxIter',1000,...
    'Display','off','Replicates',10);
% obtain the size of each cluster and the event index for each cluster
clust_size = zeros(1,cnum);
clust_ind = cell(1,cnum);
clust_ind_all = cell(1,cnum);
for c = 1:cnum
    ind_temp = find(idx_best == c);
    clust_ind_all{c} = ind_temp;
    base_up = 0;
    base_down = 0;
    base_down = length(good_ind)+base_down;
    inds = ind_temp(ind_temp > base_up & ind_temp <= base_down);
    clust_size(1,c) = length(inds);
    clust_ind{1,c} = inds;
    base_up = length(good_ind)+base_up;
end

% Obtain the amplitude / delay of the centroid 
clust_ind_origin = cell(1,cnum); % The index for the chosen ones in all the trials in the concatenated sessions
ind_flat = good_ind;
for c = 1:cnum
    clust_ind_origin{c} = [clust_ind_origin{c},ind_flat(clust_ind{c})];
end

T = 150;
[temp_delays,temp_mags,temp_true]=obtain_centroid_v2(clust_ind_origin, feat_all, T,env_true);

% sort the cluster order based on their activation map
[temp_mags_sort,temp_delays_sort,temp_true_sort,clust_ind_sort,...
    clust_ind_origin_sort,clust_size_sort] = order_clusters(clust_ind_origin,clust_ind_all,clust_size,temp_delays,temp_mags,temp_true);

% sanity check for the features
feature_ordered = feature_rescale;
base = 0;
for c = 1:cnum
    ind_temp = clust_ind_sort{c};
    feature_ordered(base + (1:sum(clust_size_sort(:,c))), :) = feature_rescale(ind_temp,:);
    base = base + length(ind_temp);
end
figure; imagesc(feature_ordered);
colorbar; xlabel('Features (16 delay + 16 power)');
ylabel('Ripple events'); set(gca,'FontSize',12,'FontWeight','bold');

%% plot the clustering result after sorting
plot_cluster_centroid(temp_mags_sort,temp_delays_sort,temp_true_sort,cnum,pos);

%% Save the results
save(['clustered_SWRs.mat'],'temp_mags_sort','temp_delays_sort','temp_true_sort','clust_ind_sort',...
    'clust_ind_origin_sort','clust_size_sort','env_true','cnum','pos','feat_all','mag_true');

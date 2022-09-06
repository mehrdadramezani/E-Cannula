function[temp_mags_sort,temp_delays_sort,temp_true_sort,clust_ind_sort,clust_ind_origin_sort,clust_size_sort]...
    = order_clusters(clust_ind_origin,clust_ind,clust_size,temp_delays,temp_mags,temp_true)
n_clust = length(clust_ind_origin);
adj_mat = zeros(n_clust,n_clust);
for n = 1:n_clust
    for m = 1:n_clust
        adj_mat(n,m) = corr2(squeeze(temp_mags(n,:,:)), squeeze(temp_mags(m,:,:)));
    end
end

adj_mat = diag(sum(adj_mat,2)) - adj_mat;
[V, D] = eig(adj_mat);

[coeff, orders] = sort(V(:,2));
temp_delays_sort = temp_delays(orders,:,:);
temp_mags_sort = temp_mags(orders,:,:);
temp_true_sort = temp_true(orders,:,:);
clust_ind_sort = clust_ind(orders);
clust_ind_origin_sort = clust_ind_origin(orders);
clust_size_sort = clust_size(:,orders);
end
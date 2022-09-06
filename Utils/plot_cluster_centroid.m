function plot_cluster_centroid(temp_mags,temp_delays,temp_true,cnum,pos)
nstep = size(temp_true,3);
figure('pos',[10+20,10,cnum*100,200]);
for k = 1:cnum
    subaxis(2,cnum,k,'SpacingVert',0.03,'MR',0.01,'SpacingHoriz',0.01,'ML',0.01,'MT',0.07,'MB',0.01);
    data = squeeze(temp_mags(k,:,:)); 
%     mag_range = [0,700];
%     mag_range = [0,2];
    mag_range = [0,0.8];
    % interpolate the data
    [X,Y] = meshgrid(1:4,1:4);
    [Xq,Yq] = meshgrid(1:0.2:4,1:0.2:4);
    dataq = interp2(X,Y,data,Xq,Yq);
    imagesc(1:0.2:4,1:0.2:4,dataq','AlphaData',~isnan(dataq')); axis square; 
    caxis(mag_range); colormap jet; % h=colorbar; caxis(mag_range); set(h,'FontSize',12); 
    title(['Cluster ',num2str(k)]); set(gca,'YTick',[],'XTick',[]);
    hold on;
    xoffsets = repmat(linspace(1,4,4),1,4);
    yoffsets = reshape(repmat(linspace(1,4,4)',1,4)',1,[])+0.2;
    t_plot = ((1:nstep)-round(nstep/2))/nstep/1.5;
    interval = 1:nstep;
    data_seg = squeeze(temp_true(k,:,:));
    data_remap = nan(16,nstep);
    data_remap(pos,:) = data_seg;
    for c = 1:size(data_remap,1)
        plot(t_plot(interval)+xoffsets(c),-data_remap(c,interval)/100+yoffsets(c),'k'); hold on;
        plot(xoffsets(c)+[0,0], [-0.4,0.3]+yoffsets(c), '--k');
    end
%     hold on;
%     data_seg = squeeze(temp_recon(k,:,:));
%     data_remap = nan(16,nstep);
%     data_remap(pos(good_ch),:) = data_seg;
%     for c = 1:size(data_remap,1)
%         plot(t_plot(interval)+xoffsets(c),-data_remap(c,interval)+yoffsets(c),'r'); hold on;
%     end
    xlim([0.7,4.3]); ylim([0.7,4.3]);
    set(gca,'XColor','None','YColor','None'); box off;
    
    subaxis(2,cnum,k+cnum,'SpacingVert',0.03,'MR',0.01,'SpacingHoriz',0.01,'ML',0.01,'MT',0.07,'MB',0.01);
    data = squeeze(temp_delays(k,:,:));
    data = data - min(data(:));
    % interpolate the data
    [X,Y] = meshgrid(1:4,1:4);
    [Xq,Yq] = meshgrid(1:0.2:4,1:0.2:4);
    dataq = interp2(X,Y,data,Xq,Yq);
    imagesc(1:0.2:4,1:0.2:4,dataq','AlphaData',~isnan(dataq')); axis square; %h=colorbar;  set(h,'FontSize',12);
    set(gca,'YTick',[],'XTick',[]); caxis([0,15]); % uo to 10 ms
    [FX,FY] = gradient(dataq');
    FX = FX/5; FY = FY/5;
    hold on; quiver(Xq(2:3:end-1,2:3:end-1),Yq(2:3:end-1,2:3:end-1),FX(2:3:end-1,2:3:end-1),FY(2:3:end-1,2:3:end-1),0,'w');
    xlim([0.7,4.3]); ylim([0.7,4.3]);
    set(gca,'XColor','None','YColor','None'); box off;
end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for "E-Cannula reveals anatomical diversity
% in sharp-wave ripples as a driver for the recruitment of distinct
% hippocampal assemblies" published in Cell Reports.
% (C) Mehrdad Ramezani, Kuzum Lab, University of California San Diego
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements threshold based ripple detection algorithm on example LFP data.
% The result (ripple time, ripple duration/amplitude/frequency) is saved to
% mat file.


% Load the data (LFP)
load('lfp_data_example.mat');

n_channels = length(lfp_data); % NEED TO CHECK NAMES (lfp_data)

% bandpass filter the data at ripple range
d_band = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',120,'HalfPowerFrequency2',250, ...
        'SampleRate',1e3);

% compute the envelope of the ripple band ECoG
ripple_data = filtfilt(d_band,data_low);
ripple_env = abs(hilbert(ripple_data));

% ripple detection parameters
thresh_dur = 20; % minimum duration is 20 ms
mask_out_len = 100; % mask out first and last 100 ms data, where ripple filtering can have transient artifacts
merge_criteria = 10; % merge adjacent ripple events within 10 ms
folds = 3.5; % multiple of SD for ripple detection
thresh_base = mean(ripple_env); % mean activity
thresh = folds*std(ripple_env,1) + thresh_base; % the threshold used to detect ripple

% detect all ripple candidates for all channels
seg_cross = ripple_env > thresh;
seg_cross(1:mask_out_len,:) = 0;
seg_cross(end-mask_out_len:end,:) = 0;


%% detect ripple events in each channel
ripple_result = struct;
for channel = 1:n_channels
    % locate the cross threshold position
    cross_start = find(diff(seg_cross(:,channel)) == 1)+1; % identify the positive cross points
    cross_end = find(diff(seg_cross(:,channel)) == -1)+1; % identify the negative cross points
    
    % Step 1 refine the start and end time of each event to baseline threshold
    for e = 1:length(cross_end)
        init_t = cross_start(e);
        end_t = cross_end(e);
        while(ripple_env(init_t,channel) > thresh_base(channel))
            init_t = init_t - 1;
        end
        while(ripple_env(end_t,channel) > thresh_base(channel))
            end_t = end_t + 1;
        end
        cross_start(e) = init_t;
        cross_end(e) = end_t;
    end
    
    % Step 2 merge the close adjacent ripples into one (closer than threshold inter-duration)
    t_inter = cross_start(2:end) - cross_end(1:end-1);
    ind = t_inter < merge_criteria;
    cross_start(find(ind)+1) = [];
    cross_end(ind) = [];
    
    % Step 3 throw away the events with duration < threshold
    dur = cross_end - cross_start;
    ind_keep = dur >= thresh_dur;
    cross_start = cross_start(ind_keep);
    cross_end = cross_end(ind_keep);
    
    ripple_result(channel).ripple_start = cross_start;
    ripple_result(channel).ripple_end = cross_end;
    ripple_result(channel).ripple_dur = cross_end - cross_start;
    
end

% Now obtain the ripple event across array by merging all channels
combine_ripple = t_low.*0;
for channel = 1:n_channels
    for i = 1:length(ripple_result(channel).ripple_start)
        start = ripple_result(channel).ripple_start(i);
        endt  = ripple_result(channel).ripple_end(i);
        combine_ripple(start:endt)=1;
    end
end
ripple_result_combine = struct;
ripple_result_combine.ripple_start = find(diff(combine_ripple)==1)'+1; % sample index at 1 kHz
ripple_result_combine.ripple_end = find(diff(combine_ripple)==-1)';
ripple_result_combine.ripple_dur = ripple_result_combine.ripple_end - ...
    ripple_result_combine.ripple_start;
ripple_result_combine.ripple_center = ceil(0.5*(ripple_result_combine.ripple_end + ...
    ripple_result_combine.ripple_start));

% Now start manual curation for the detected ripple events one-by-one...
% now go through all the ripple events and label them...
% Instructions: 
% - Use left and right keys to go between different ripple events
% - Use up key to mark it true ripple
% - Use down key to mark it as false ripple

lfp_mappped(:,pos) = lfp_data;

true_ripples = [];
current = 1;
fig1 = figure('pos',[700,10,300,700]);
hax1 = axes; hold on;
fig3 = figure('pos',[1000,110,600,600]);
hax3 = axes; hold on;
fig2 = figure('pos',[100,110,600,600]);
hax2 = axes;
while(current >= 1 && current <= n_all)
%     if ismember(current, throw)
%         current = current + 1;
%         continue;
%     end
    
    inds = (ripple_result_combine.ripple_start(current)-100) : (ripple_result_combine.ripple_end(current)+100);
    if min(inds)<1
        current = current + 1;
        continue;
    end
    data = ripple_env(inds,:);
    
    data_amp = sqrt(mean(data(101:end-100,:).^2));
    data_amp(pos) = data_amp;
    data_amp = reshape(data_amp,4,4);
    mag_range = [0,prctile(data_amp(:),90)];
    % interpolate the data
    [X,Y] = meshgrid(1:4,1:4);
    [Xq,Yq] = meshgrid(1:0.2:4,1:0.2:4);
    dataq = interp2(X,Y,data_amp,Xq,Yq);
    imagesc(hax2,1:0.2:4,1:0.2:4,dataq','AlphaData',~isnan(dataq')); axis square; h=colorbar; caxis(mag_range);
    set(h,'FontSize',12); set(gca,'YTick',[],'XTick',[]); 
    hold on;
    xoffsets = repmat(linspace(1,4,4),1,4);
    yoffsets = reshape(repmat(linspace(1,4,4)',1,4)',1,[])+0.2;
    nstep = length(inds);
    t_plot = ((1:nstep)-round(nstep/2))/nstep/1.2;
    interval = 1:nstep;
    data_remap = nan(16,nstep);
    data_remap(pos,:) = data';
    data_remap = data_remap /100; % rescale to fit into the plot
    for c = 1:size(data_remap,1)
        plot(hax2,t_plot(interval)+xoffsets(c),-data_remap(c,interval)+yoffsets(c),'k'); hold on;
        plot(hax2,xoffsets(c)+[0,0], [-0.4,0.3]+yoffsets(c), '--k');
        plot(hax2,t_plot(interval(101:end-100))+xoffsets(c), -data_remap(c,interval(101:end-100))+yoffsets(c),'r');
    end
    
    [~,ind] = min(abs(t_Ca - ripple_result_combine.ripple_start(current)/fs_low));
    if move_ind(ind)
        title(hax2,['Ripple ',num2str(current),'/',num2str(n_all),' move']);
    else
        title(hax2,['Ripple ',num2str(current),'/',num2str(n_all)]);
    end
    xlim(hax2,[0.5,4.5]); ylim([0.5,4.5]);
    set(hax2,'XColor','None','YColor','None'); box off;
    hold off;
    
    % plot the LFP around it
    ind_plot = ripple_result_combine.ripple_center(current);
    for i = 1:16
        plot(hax1,t_low(ind_plot+(-100:100)), lfp_mappped(ind_plot+(-100:100),i)+i*500); hold on;
        text(hax1,t_low(ind_plot+80), i*500+50, num2str(i),'fontsize',12);
    end
    plot(hax1,[1,1]*t_low(ind_plot),[500,8000],'r');
    hold off;
    xlim(hax1,[-0.1,0.1]+ind_plot/fs_low);
    ylim(hax1,[0,9000]);
    
    
    
    % Plot the LFP around it (mapped)
    t_plot = t_low(ind_plot+(-100:100)) - t_low(ind_plot);
    for i = 1:16
        plot(hax3,t_plot+xoffsets(i)/3,-lfp_mappped(ind_plot+(-100:100),i)/4000+yoffsets(i)/3,'color',[0, 0.4470, 0.7410]); hold on;
        text(hax3,t_plot(1)+xoffsets(i)/3-0.05, yoffsets(i)/3, num2str(i),'fontsize',12);
        plot(hax3,xoffsets(i)/3+[0,0], [-0.1,0.1]+yoffsets(i)/3, '--r');
    end
    %hold off;
    xlim(hax3,[0.1,1.5]); ylim(hax3, [0.1,1.5]);
    set(hax3,'XColor','None','YColor','None'); box off;
    set(hax3,'YDir','reverse');
    
    k = waitforbuttonpress; 
    value = double(get(fig2,'CurrentCharacter'));
    if value == 30 % up
        true_ripples = [true_ripples, current];
        true_ripples = unique(true_ripples);
        disp(['Added ripple ', num2str(current)]);
    elseif value == 31 % down
        true_ripples = setdiff(true_ripples, current);
        disp(['Removed ripple ', num2str(current)]);
    elseif value == 28 % left
        current = current - 1;
    elseif value == 29 % right
        current = current + 1;
    else
        break;
    end
    cla(hax1);
    cla(hax2);
    cla(hax3);
end


%  Save the ripple detection results
save(['refined_SWRs.mat'],'true_ripples','ripple_result_combine');

function [detected_spike_locs, spikesBelowThresh, spikeTemplate, spikeDistThreshold, filter_params] = extracellular_spike_extractor_03_02_2016(...
    unfiltered_data, piezo, template, spikeTemplateWidth,  numtemps, Distance_threshold, approval, approval_time)
%% Code to match a given spike waveform across a recording (1-D); basic idea is to estimate warping distance via DTW and convert the  estimated warp distance into a probability via a kernel
last_filt_params = 'C:\Data\spikes\last_filt_params.mat';
if exist(last_filt_params, 'file');load(last_filt_params);end

global vars;
vars.fs = 10000;
vars.spikeTemplateWidth = spikeTemplateWidth;

spikeDistThreshold  = Distance_threshold;

max_len = 400000;
if length(unfiltered_data) < max_len 
vars.len = length(unfiltered_data)-1000;
else
vars.len = max_len -1000;
end

start_point = find(unfiltered_data >  mean(unfiltered_data)+ 3*std(unfiltered_data),50000);
stop_point = start_point+max_len;

%% choose the filter parameters for pulling out spikes
filter_sliderGUI(unfiltered_data(start_point:stop_point), piezo(start_point:stop_point),spikeTemplateWidth);

if ~isempty(template)
spikeTemplate = template;
spikeTemplateWidth = length(spikeTemplate);
else
    
fig = gcf;dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')

spikes_flag = [];
while isempty(spikes_flag)          % Wait while the user does this.
   
    reuse = input('\n To use last spike template and HWD, press 1, otherwise hit enter...');
if reuse ==1 
       spikeTemplate = vars.spikeTemplate;
       Distance_threshold = vars.Distance_threshold;
else
  for hh = 1:numtemps
        spikes_flag = input('\n click on a spike peak and hit enter...');
    try
        fig;template_center = [];
        c_info = getCursorInfo(dcm_obj);
        template_center = c_info.Position(1);
        spikeTemplate(hh,:) = vars.filtered_data(round(template_center-(0.5*spikeTemplateWidth)):round(template_center+(0.5*spikeTemplateWidth)));
        catch err
        display('no spikes selected'); detected_spike_locs = []; spikesBelowThresh = []; spikeTemplate = template;
    end
   end
end
    spikes_flag = 0;
        if size(spikeTemplate,1)>1; spikeTemplate = mean(spikeTemplate,1);end
end

if isempty(spikeTemplate); detected_spike_locs = []; spikesBelowThresh = []; 
return; end;

close;
figure(12), plot(spikeTemplate), title('Template Waveform: wait one second');pause(1);close;
vars.spikeTemplate = spikeTemplate;
end

%% get all the spike locs using the correct filt and thresh cvalues
spikeTemplateWidth = length(spikeTemplate);
filts = vars.hp_cutoff/(vars.fs/2);
[x,y] = butter(4,filts,'high');%%bandpass filter between 50 and 200 Hz
filtered_data_high = filter(x, y, unfiltered_data);

filts = vars.lp_cutoff/(vars.fs/2);
[x2,y2] = butter(4,filts,'low');%%bandpass filter between 50 and 200 Hz
filtered_data = filter(x2, y2, filtered_data_high);

if vars.diff == 0
diff_filt = filtered_data';
elseif vars.diff == 1
diff_filt = [0 diff(filtered_data)'];
diff_filt(1:100) = 0;
elseif vars.diff == 2
diff_filt = [0 0 diff(diff(filtered_data))'];
diff_filt(1:100) = 0;
end

all_filtered_data = diff_filt;

[locks, ~] = peakfinder(double(all_filtered_data),mean(all_filtered_data)+vars.peak_threshold*std(all_filtered_data));%% slightly different algorithm;  [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local maxima that are at least sel above surrounding all_filtered_data and larger (smaller) than thresh if you are finding maxima (minima).
loccs = locks((locks> spikeTemplateWidth)); %%to prevent strange happenings, make sure that spikes do not occur right at the edges
spike_locs = loccs(loccs < (length(all_filtered_data)-spikeTemplateWidth));


%% pool the detected spike candidates and do template matching
targetSpikeDist = zeros(size(spike_locs));
counter = 1;
for i=1:length(spike_locs)

    if min(spike_locs(i)+spikeTemplateWidth/2,length(all_filtered_data)) - max(spike_locs(i)-spikeTemplateWidth/2,0)< spikeTemplateWidth
        continue
    else
%         curSpikeTarget = all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2),0)+1:min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)),1);  %% old method
%         detectedSpikeCandidates(:,counter) = all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2),0): min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)),1);
        curSpikeTarget = all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2),0)+1: min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)));
        detectedSpikeCandidates(:,counter) = all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2),0): min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)));
        norm_curSpikeTarget = curSpikeTarget/max(curSpikeTarget);norm_spikeTemplate = spikeTemplate/max(spikeTemplate);
        % figure(34); hold all; plot(curSpikeTarget, 'b');plot(norm_curSpikeTarget, 'r');pause;
        [targetSpikeDist(i), foo,bar] = dtw_WarpingDistance(norm_curSpikeTarget, norm_spikeTemplate);
        clear foo bar;      
    end
    counter = counter+1;
end

% estimate spike probabilities at candidate locations
spikeProbs = zeros(size(spike_locs));
for i=1:length(spike_locs)
   if min(spike_locs(i)+spikeTemplateWidth/2,length(all_filtered_data)) - max(spike_locs(i)-spikeTemplateWidth/2,0)< spikeTemplateWidth
       continue
   else
       spikeProbs(i) = exp( -(abs(targetSpikeDist(i)-mean(targetSpikeDist))) / (2*var(targetSpikeDist)) );
   end
end

%%
figure(11); clf; set(11, 'Position', [0 0 1000 700],'color', 'w');
accept_spikeDistThreshold_flag = 0;
while accept_spikeDistThreshold_flag == 0
detected_spike_locs = [];clear spikesAboveThresh spikesBelowThresh spikesAboveThresh_aligned spikesBelowThresh_aligned;
figure(11); clf;
subplot(3,2,1); plot(detectedSpikeCandidates), title('All detected spike candidates');hold on; plot(spikeTemplate, 'k', 'linewidth', 4)
subplot(3,2,3), hist(targetSpikeDist,20), title('Histogram of distances')
subplot(3,2,5);plot(unfiltered_data-mean(unfiltered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro');axis off;

    spikesAboveThresh = []; spikesBelowThresh = [];spikeDistThreshold = [];detected_spike_locs = [];
    
    if isempty(Distance_threshold)
    subplot(3,2,5);title('Enter the new spike distance threshold:');hold all;
    spikeDistThreshold = input('\n \n Enter the spike distance threshold:');
    vars.Distance_threshold = spikeDistThreshold;
    else
    spikeDistThreshold  = Distance_threshold;
    vars.Distance_threshold = Distance_threshold;
    end
    i =1;
    aboveThreshCounter =1;
    belowThreshCounter=1;
    
    while i<= length(spike_locs)
                if length(all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2))+1: min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)))) >= spikeTemplateWidth-1
            curSpikeTarget = all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2))+1: min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)));
                    %         if length(all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2),0)+1: min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)),1)) >= spikeTemplateWidth-1
%             curSpikeTarget = all_filtered_data(max(spike_locs(i)-floor(spikeTemplateWidth/2),0)+1: min(spike_locs(i)+floor(spikeTemplateWidth/2),length(all_filtered_data)),1);
            if targetSpikeDist(i)> spikeDistThreshold
                spikesAboveThresh(:, aboveThreshCounter)= curSpikeTarget;
                aboveThreshCounter = aboveThreshCounter+1;
            else
                spikesBelowThresh(:, belowThreshCounter)= curSpikeTarget;
                belowThreshCounter = belowThreshCounter+1;
                detected_spike_locs = [detected_spike_locs;spike_locs(i)];
            end
        end
        i=i+1;
    end
    % update the counters
    belowThreshCounter = belowThreshCounter-1;
    aboveThreshCounter = aboveThreshCounter-1;
    
%     ind = detected_spikes<spikeDistThreshold;
%     finalSpikeLocs = detected_spike_locs;
    
    figure(11);
    subplot(3,2,2), plot(spikesAboveThresh), title(sprintf('Spikes above %f: number = %d', spikeDistThreshold, size(spikesAboveThresh,2)))
    subplot(3,2,4), plot(spikesBelowThresh), title(sprintf('Spikes below %f: number = %d', spikeDistThreshold, size(spikesBelowThresh,2)))
    subplot(3,2,6), plot(all_filtered_data-mean(all_filtered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro'); title ('Estimated spikes');
    if approval== 0
      accept_spikeDistThreshold_flag = 1;
    elseif approval == 1
       subplot(3,2,5), axis off; title({'Are you happy with the spike output?' 'Hit 0 if No'}); hold all;
       plot(unfiltered_data-mean(unfiltered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro');
       if ~exist('approval_time','var')
           approval_time = 1;
       end
       accept_spikeDistThreshold_flag = timeinput(approval_time); %%returns 1 if no input within 1 second
       Distance_threshold = [];
    elseif approval == 2
    subplot(3,2,5), axis off; title({'Are you happy with the spike output?' 'Hit Enter for yes or press 1 for No'}); hold all;
    plot(unfiltered_data-mean(unfiltered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro');
    accept_spikeDistThreshold_flag = input('\n Are you happy with the spike output? \n Hit Enter for yes or press 0 for No:');
    Distance_threshold = [];
    end
end

filter_params = vars;
save(last_filt_params, 'vars');display('spike filter parameters saved');
end

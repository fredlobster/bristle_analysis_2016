
clear all;

%% init file and paths
% files = {'2016_03_02_0006.abf', '2016_03_02_0004.abf', '2016_03_02_0005.abf', '2016_03_02_0009.abf'};
% ranges = {[1:3.8e6],[],[],[],[],[],[],[]};

files = {'2016_02_26_0008.abf', '2016_02_26_0007.abf', '2016_02_26_0006.abf', '2016_02_26_0005.abf'};
ranges = {[],[],[],[],[],[],[],[]};
gaps = [20, 40, 80, 160];
output_figure = 0;
extra = 1000;

rootfolder = 'C:\Data\';
datapath = [rootfolder 'ephys_data\'];
spikepath = [rootfolder 'spikes\'];
stimpath = [rootfolder 'stim_metadata\'];
summpath = [rootfolder 'ephys_metadata\'];
outputpath = [rootfolder 'figures\'];

fig2 = figure(2);clf;hold all;set(fig2,'color', 'w', 'position', [1 1 800 1050]);hold all;
cmap = colormap(winter(length(files)));

for abf = 1: length(files) %% if averaging is to be done across files  
filename = files{abf};clear Data;
metadata = readtable([stimpath files{1}(1:end-4) '.txt'],'Delimiter','\t');
num_conditions = height(metadata);
num_cycles = mode(table2array(metadata(1:num_conditions,7)));
num_reps = mode(table2array(metadata(1:num_conditions,5)));
expected_voltages = table2array(metadata(1:num_conditions,8));

load([summpath filename(1:end-4) '_summ']);

%% plot data
subplots = reshape(1:num_conditions*num_cycles,num_conditions, num_cycles);

ii = 1;
for cond=1:num_conditions
    for cycle = 1:num_cycles
        for rep = 1:size(Data(cond,cycle).bristle_spike_rate,1)
        subplot(num_conditions, num_cycles, ii);box off;hold on;
%         plot(Data(cond,cycle).bristle_spike_rate(rep,:), 'color', cmap(abf,:));
%         ylim([0 275]);

        plot(Data(cond,cycle).bristle_1(rep,:)-mean(Data(cond,cycle).bristle_1(rep,extra-100:extra)), 'color', cmap(abf,:));
        ylim([-25 0]);
        
        peak_firing_rates(cond,cycle,rep) = max(Data(cond,cycle).bristle_spike_rate(rep,:));
        peak_vm(cond,cycle,rep) = min(Data(cond,cycle).bristle_1(rep,extra:end))-mean(Data(cond,cycle).bristle_1(rep,extra-100:extra));
        spike_counts(cond,cycle,rep) = sum(Data(cond, cycle).bristle_spikes_1(rep,:));
%         plot(Data(abf,cond,cycle).patch_1(rep,:)-mean(Data(abf,cond,cycle).patch_1(rep,:)), 'color', cmap(rep,:));box off;hold on;
%         plot(Data(abf,cond,cycle).piezo_1(rep,:), 'color', cmap(rep,:));box off;hold on;ylim([0 10]);
        end
        ii = ii +1;
    end
end

mean_rates(abf,:,:) = sum(peak_firing_rates,3) ./ sum(peak_firing_rates~=0,3);
mean_counts(abf,:,:) = sum(spike_counts,3) ./ sum(spike_counts~=0,3);
mean_vm(abf,:,:) = sum(peak_vm,3) ./ sum(peak_vm~=0,3);


end


fig3 = figure(3);clf;hold all;set(fig3,'color', 'w', 'position', [1 1 1000 1200]);hold all;
for cond=1:num_conditions
for abf = 1:length(files)
subplot(3,6,cond);box off; hold all;
plot(squeeze(mean_counts(abf,cond,:)),'-o', 'markersize', 3, 'color',cmap(abf,:)); 
if cond ==1; ylabel('spike count');xlabel('cycle #');else set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); end
ylim([min(mean_counts(:))-5 max(mean_counts(:))+5]);
title({ [num2str(table2array(metadata(cond,2))) ' ms ON']; [num2str(table2array(metadata(cond,3))) ' ms OFF']});

subplot(3,6,cond+6);box off; hold all;
plot(squeeze(mean_rates(abf,cond,:)),'-o','markersize', 3, 'color',cmap(abf,:)); 
if cond ==1; ylabel('peak firing rate');xlabel('cycle #');else set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); end
ylim([min(mean_rates(:))-10 max(mean_rates(:))+10]);
title({ [num2str(table2array(metadata(cond,2))) ' ms ON']; [num2str(table2array(metadata(cond,3))) ' ms OFF']});

subplot(3,6,cond+12);box off; hold all;
plot(squeeze(mean_vm(abf,cond,:)),'-o','markersize', 3, 'color',cmap(abf,:)); 
if cond ==1; ylabel('peak firing rate');xlabel('cycle #');else set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); end
ylim([min(mean_vm(:))-3 max(peak_vm(:))+3]);
title({ [num2str(table2array(metadata(cond,2))) ' ms ON']; [num2str(table2array(metadata(cond,3))) ' ms OFF']});


end
end

legend('20 s', '40 s', '80 s', '160 s', 'location', 'best');

if output_figure == 1
    try export_fig(fig2,[outputpath '/' filename(1:end-2) '.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    export_fig(fig4,[outputpath '/' filename(1:end-2) '_raster.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    export_fig(fig5,[outputpath '/' filename(1:end-2) '_means.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');

    catch;end;
end
% datestr(now, 'mm_dd_yy')
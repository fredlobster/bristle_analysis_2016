CA; close all;
% day = datestr(now, 'mm_dd_yy');
day = '14_03_21';

%%pull out most recent data file
fid = fopen(['C:\LabviewRoot\Metadata\metadata_' day '.txt']);

strs = textscan(fid, '%s');
idx = strfind(strs{1}, 'data');
gal = strfind(strs{1}, 'gal4');
trues = find(~cellfun(@isempty,idx));
truesg = find(~cellfun(@isempty,gal));

for jj = 1:length(trues)
filename = strs{1}{trues(jj)};
pos_filename = ['pos' filename(5:end-2) '.txt'];
try gal4 = strs{1}{truesg(jj)};
catch err; gal4 = 'random';
    %%
end

% filename = 'data_27_09_12_150807.o';
crange = {[]};brist = 0;freqs = 1;

fileroot = 'C:\LabviewRoot\Data_files\';
fs = 10000;
data = ReadAllDataFromFile([fileroot filename]);
patch_1 = data(2,:)./3; %%patch primary signal
patch_2 = data(4,:)./3*100; %%patch current signal
bristle_1 = data(1,:)./3; %%bristle primary signal
bristle_1 = bristle_1-mean(bristle_1)-27;
bristle_2 = data(3,:)./3*100; %%bristle current signal (1 mv/ pA)
bristle_2 = (bristle_2-mean(bristle_2))/8-50;
% patch_3 = data(6,:)./3*100; %%bristle command potential


piezo_1 = -data(5,:)./4;
piezo_1 = (piezo_1-mean(piezo_1));
piezo_2 = data(6,:)./4;
piezo_2 = (piezo_2-mean(piezo_2));
patch_1_range = max(patch_1) - min(patch_1);

time_vec = (1:length(patch_1))/fs;
secs = length(patch_1)/fs;
plot_range = 1:secs*fs;

txtfile = [fileroot '\' pos_filename];
txt_filename = txtfile;
try
values = read_stepper_text_12_2013(txtfile);
x_pos = values(2,:)-mean(values(2,:));x_pos = x_pos-mean(x_pos);
y_pos = values(3,:)-mean(values(3,:));y_pos = y_pos-mean(y_pos);
z_pos = values(4,:)-mean(values(4,:)); z_pos = z_pos-mean(z_pos);
x_pos2 = values(5,:)-mean(values(5,:));x_pos2 = x_pos2-mean(x_pos);
y_pos2 = values(6,:)-mean(values(6,:));y_pos2 = y_pos2-mean(y_pos);
z_pos2 = values(7,:)-mean(values(7,:)); z_pos2 = z_pos2-mean(z_pos);
stepper_failed = [];
catch stepper_failed;
end
 
stepper_scaling = 50;
offset = -max(patch_1)+2;
piezo_offset = -min(patch_1)+2;
piezo_scaling = 6;

if range(piezo_1) > .25
    piezo  = piezo_1;
elseif range(piezo_2) > .25
    piezo  = piezo_2;
else 
    piezo = zeros(1,length(patch_1));
end

fig1 = figure(1);clf;hold all;set(fig1,'color', 'w', 'position', [1 1 800 1050]);
subplot(3,1,1); hold all;box off;
plot(time_vec(plot_range),patch_1(plot_range), 'k'); 
% plot(time_vec(plot_crange),filtered_data(plot_crange)+ max(patch)+1, 'color' , [0.5 0.5 0.5]); 
plot(time_vec(plot_range), piezo_1(plot_range)/piezo_scaling-piezo_offset, 'r');
plot(time_vec(plot_range), piezo_2(plot_range)/piezo_scaling-piezo_offset-1, 'b');
if isempty(stepper_failed)
plot(time_vec((1:1000:end)),y_pos(1:ceil(length(time_vec)/1000))*stepper_scaling-offset, 'm');
plot(time_vec((1:1000:end)),y_pos2(1:ceil(length(time_vec)/1000))*stepper_scaling-offset, 'g');
end
% if secs > 31;xlim([0 30]);end
title({[gal4 ', ' filename]}, 'interpreter', 'none');
axis tight;

try
decs = find(diff(piezo) < -0.25);
incs = find(diff(piezo) > 0.25);
if decs(1) < incs(1);decs(1) = [];end
if any(diff(incs) < 100); incs(diff(incs) < 10) = [];end %%make sure incs are not too close
if any(diff(decs) < 100); decs(diff(decs) < 10) = [];end
if incs(length(incs)) > decs(length(decs)); incs(length(incs)) = [];end

per = 2.*abs(decs(1)-incs(1));
extra = per/10;
patches_1 = [];stims = [];ii = 1;patches_2 = [];bristles_1 = [];bristles_2 = [];
for ii = 1:length(decs)
if incs(ii)+per < length(patch_1) && incs(ii)-extra > 0
    patches_1(ii,:) = patch_1(incs(ii)-extra:incs(ii)+per)';
    patches_2(ii,:) = patch_2(incs(ii)-extra:incs(ii)+per)';
    stims(ii,:) = piezo(incs(ii)-extra:incs(ii)+per)';
    bristles_1(ii,:) = bristle_1(incs(ii)-extra:incs(ii)+per)';
    bristles_2(ii,:) = bristle_2(incs(ii)-extra:incs(ii)+per)';
ii = ii+1;
end
end

cmap = colormap(autumn(length(decs)));
subplot(3,1,2); hold all;box off;

for jj = 1:length(decs)
ext = time_vec(1:decs(jj)-incs(jj)+401);
plot(ext, patch_1((incs(jj)-200):(decs(jj)+200)), 'linewidth', 0.75, 'color', cmap(jj,:))
plot(ext, piezo((incs(jj)-200):(decs(jj)+200))/piezo_scaling-piezo_offset, 'b', 'linewidth', 0.75);
end
axis tight;

subplot(3,1,3);hold all;box off;
confplot((1:size(patches_1,2))/fs,mean(patches_1), std(patches_1), std(patches_1), 'k');hold on;
if brist ==1; confplot((1:size(bristles_2,2))/fs,mean(bristles_2), std(bristles_2), std(bristles_2), 'b');hold on;end
legs(1) = plot((1:size(patches_1,2))/fs,mean(stims)/piezo_scaling-piezo_offset, 'b', 'linewidth', 0.75);hold on;

% legs(2) = plot((1:size(patches_1,2))/fs,mean(bristles_1), 'b', 'linewidth', 0.75);
% legs(3) = plot((1:size(patches_1,2))/fs,mean(patches_1), 'k', 'linewidth', 0.75);
% legend(legs,'piezo signal', 'bristle signal', 'V_m', 'location', 'northeast');

% ylim([-60 0]);
ylabel('V_m (mV)');
xlabel('time (s)');
axis tight;
% h =  text(0.52,5, );set(h, 'FontWeight', 'bold', 'Fontsize', 16) 
catch err
end

file_root = 'C:\MatlabRoot\all_figures\';
if ~exist([file_root day],'dir')
mkdir([file_root day ]);
end
[SUCCESS,MESSAGE,MESSAGEID]= mkdir([file_root]);
% export_fig(fig2,['C:\MatlabRoot\central_recordings\10_24_2012' '/' file_name '_averages.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
export_fig(fig1,[file_root day '\' filename(1:end-2) '_figure.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
end
% % % % 
close all;
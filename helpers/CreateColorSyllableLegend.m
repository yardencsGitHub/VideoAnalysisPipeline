%% create color legend and syllable spectrograms
separation = 0.05; % seconds, image separation
freq_min = 00; freq_max = 8000;

bird_name = 'lbr3022'; %'lrb85315';
bird_folder_name = 'lbr3022'; %'lrb853_15';
template_file = 'lbr3022_template'; %'lrb85315template';
%annotation_file = 'lrb85315auto_annotation5';
CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

syllables = [5 6 7]; %[templates.wavs.segType]; %[0 1 2 8 207 300 301 304]; %[0:4 8 9 200:209 300 301 303 304]; 
%syllables = syllables((idx(22:26)));
%l[5 6 7 302]
% [0 1 2 8 207 300 301 304]
%setdiff([0:9 200:209 300:309 400:409 500 ],[400 500 2 200 202 204 300 301 302 303 0 1 4 5 8 9 203 208 209 304:306 308:309 408]);
%setdiff([0:9 200:209 300:309 400:409 500 ], [0 1 4 5 8 9 203 208 209 300:306 308:309 408]); %[0 1 2 4 5 8 9 200 202:204 208 209 300:306 308:309 408];%[0:9 200:209 300:309 400:409 500];
extras = 3;
%[400 409 500 2 200 202 204 300 301 302 303]
%[0:2 4 5 8 9 200 202:204 208 209 300:306 308:309 400 408 500]
%[0 1 4 5 8 9 203 208 209 304:306 308:309 408]

%%
dt = 0.9167*0.001;
FS = 48000;
fbins = 257;
separation = floor(separation/dt)*dt;

addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/zftftb'));

cd(laptop_wav_folder);
load(fullfile(laptop_wav_folder,[template_file '.mat']));
tmp_syllables = [[templates.wavs.segType] ]; %-1 102 103
n_syllables = numel(syllables);

colors = distinguishable_colors(numel(tmp_syllables));
n_bins = 0;
IMAGE = [];
IMAGE1 = [];
h = figure;
subplot(21,1,1);
curr_start = 0;
sig_lns = [];
for sylnum = 1:numel(templates.wavs)
    if ismember(templates.wavs(sylnum).segType,syllables)
        n_bins = n_bins + ceil(numel(templates.wavs(sylnum).wav)/FS/dt+separation/dt);
        sig = templates.wavs(sylnum).wav;
        
        [S,F,T,P] = spectrogram((sig/(sqrt(mean(sig.^2)))),220,220-44,512,FS,'reassigned');
        sig_len = size(S,2);
        sig_lns = [sig_lns sig_len];
        IMAGE = [IMAGE P zeros(fbins,ceil(separation/dt))]; %abs(S)
        IMAGE1 = [IMAGE1 abs(S) zeros(fbins,ceil(separation/dt))]; %abs(S)
        plot(curr_start+dt:dt:curr_start+sig_len*dt,ones(1,sig_len),...
            'LineWidth',10,'Color',colors(find(tmp_syllables == templates.wavs(sylnum).segType),:));
        text(curr_start,2,num2str(templates.wavs(sylnum).segType),'FontSize',16,'Color',[0 0 0]);
        curr_start = curr_start + sig_len*dt + separation;
        
        hold on;
    end
end
xlim([0 curr_start])
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'color','none');
axis off;
ylim([0.9 2]);
%%
subplot(21,1,2:11);
imagesc(dt:dt:curr_start,F,IMAGE1);

xlim([0 curr_start])
set(gca,'Ydir','normal');
colormap(hot)
ylim([freq_min freq_max]);    
hold on;
tmp = 0;

for cnt = 1:numel(sig_lns)
        line([tmp tmp],[freq_min freq_max],'Color',[1 1 1]/2,'LineStyle','--');
        line([tmp+sig_lns(cnt)*dt tmp+sig_lns(cnt)*dt],[freq_min freq_max],'Color',[0.5 0.5 0.5],'LineStyle','--');
        
        tmp = tmp + sig_lns(cnt)*dt + separation;
end
set(gca,'XTick',[]);

set(gca,'YTick',[])
caxis([0 100]);
%%
subplot(21,1,12:21);
imagesc(dt:dt:curr_start,F,log(IMAGE*50000+1));

xlim([0 curr_start])
set(gca,'Ydir','normal');
colormap(hot)
ylim([freq_min freq_max]);
caxis([0.5 5.5])  
tmp = 0;
xt = [];
for cnt = 1:numel(sig_lns)
        line([tmp tmp],[freq_min freq_max],'Color',[1 1 1]/2,'LineStyle','--');
        line([tmp+sig_lns(cnt)*dt tmp+sig_lns(cnt)*dt],[freq_min freq_max],'Color',[0.5 0.5 0.5],'LineStyle','--');
        xt = [xt tmp];
        tmp = tmp + sig_lns(cnt)*dt + separation;
end
xtickformat('%2.0f')
set(gca,'XTick',xt+sig_lns(cnt)*dt/2);
tlbl = {};
for cnt = 1:numel(sig_lns)
    tlbl = {tlbl{:} ['\color{black}' sprintf('%3.0f',sig_lns(cnt))]};
end
set(gca,'XTickLabel',tlbl); %sig_lns*dt*1000);


set(gca,'YTick',[500 1500 3000 4500],'Color','k')
set(gca,'FontSize',16);
xlabel('Syllable durations (msec)','Color',[0 0 0]);
ylabel('frequency (Hz)','Color',[0 0 0]);
colormap(copper)
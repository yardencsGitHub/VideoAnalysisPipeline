%%
% Folders on laptop:
laptop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs';
laptop_wav_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav';
laptop_gif_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/gif';
laptop_storage_folder = '/Volumes/CanaryData/DATA/lrb853_15/movs/';
laptop_annotated_dir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated';
laptop_annotated_images_dir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images';
DamagedFolder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/too_large_or_damaged/';
laptop_manualROI_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs';

% Folders on Data desktop:
desktop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs';
desktop_storage_folder = '/Volumes/home/Data/Imaging/lrb853_15/RawData';
desktop_max_projections_dir = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs/MaxProj';

% Parameters
last_idx = 7982;
last_date = '2017_06_29';
last_time = '08_24_46';
bird_name = 'lrb85315';

%% Create single song phrase annotated ROI (zscored) activity
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/zftftb'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
cd (laptop_manualROI_folder);
load lrb85315template;
syllables = [[templates.wavs.segType] -1 102 103];
load lrb85315auto_annotation5;
ord = [];
dates = [];
for i = 1:numel(keys)
    tokens = regexp(keys{i},'_','split');
    ord = [ord; str2num(tokens{2})];
    dates = [dates; char(join(tokens(3:5),'_'))];
end
[locs,indx] = sort(ord);
elements = elements(indx);
keys = keys(indx);
dates = dates(indx,:);
cd ROIdata;
%%
n_del_frames = 4;

for cnt = 963:numel(ord)
    cd([laptop_manualROI_folder '/ROIdata']);
    hashnum = ord(cnt); %8110;

    loc = min(find(locs == hashnum));
    phrases = return_phrase_times(elements{loc});

    cd(dates(loc,:)); % 2017_07_06;
    try
        load(['ROIdata_' keys{loc}(1:end-3) 'mat']); % lrb85315_8105_2017_07_06_05_37_56;
    catch em
        if strcmp(dates(cnt,:),'2017_04_19')
            continue;
        end
        tokens = regexp(keys{loc},'_','split');
        tokens{2} = num2str(str2num(tokens{2})); %sprintf('%04d',str2num(tokens{2}));
        fname = char(join(tokens,'_'));
        [SUCCESS,MESSAGE,MESSAGEID] = movefile(['ROIdata_' fname(1:end-3) 'mat'],['ROIdata_' keys{loc}(1:end-3) 'mat']);
        load(['ROIdata_' keys{loc}(1:end-3) 'mat']);
        display('fixed file name');
    end


    n_syllables = numel(syllables);
    freq_min = 300; freq_max = 8000;
    colors = distinguishable_colors(n_syllables);

    h=figure('Visible','off','Position',[77          91        2215         420]);

    subplot(10,1,1);
    t = vidTimes;    
    for phrasenum = 1:numel(phrases.phraseType)
        tonset = phrases.phraseFileStartTimes(phrasenum);
        toffset = phrases.phraseFileEndTimes(phrasenum);
        plot(t((t>=tonset) & (t<=toffset)),ones(1,sum((t>=tonset) & (t<=toffset))),'Color', ...
            colors(find(syllables == phrases.phraseType(phrasenum)),:),'LineWidth',10);
        hold on;
    end
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    xlim([vidTimes(n_del_frames+1) t(end)]);
    set(gca,'color','none');
    axis off;
    ylim([0.9 2]);
    tokens = regexp(keys{loc}(1:end-3),'_','split');
    title(['bird: ' tokens{1} ', file: ' tokens{2}]);
    
    subplot(10,1,2:10);
    s = zscore(dff(:,n_del_frames+1:end)')';
    imagesc(vidTimes(n_del_frames+1:end),1:size(s,1),s);
    colormap(hot); caxis([quantile(s(:),0.75) max(s(:))]);

    for phrasenum = 1:numel(phrases.phraseType)
            tonset = phrases.phraseFileStartTimes(phrasenum);
            toffset = phrases.phraseFileEndTimes(phrasenum);
            line([tonset tonset],[0 size(s,1)+1],'Color',[1 1 1],'LineStyle','--');
            line([toffset toffset],[0 size(s,1)+1],'Color',[0.5 0.5 0.5],'LineStyle','--'); 
    end
    
    
    saveas(h,fullfile(laptop_manualROI_folder,'SingleFileImages',[keys{loc}(1:end-3) 'png']));
    hgclose(h);
end
%cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs/ROIdata');

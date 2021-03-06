%%
% Parameters
%last_idx = 7982;
%last_date = '2017_06_29';
%last_time = '08_24_46';
raw_data_prefix = 'baseROIdata_'; %'NonoverlapBaseROIdata_'

bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009auto_annotation1_fix'};
bird_params = bird3_params;
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
% Folders on laptop:
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

% Folders on Data desktop:
desktop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs';
desktop_storage_folder = '/Volumes/home/Data/Imaging/lrb853_15/RawData';
desktop_max_projections_dir = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs/MaxProj';



%% Create single song phrase annotated ROI (zscored) activity
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/zftftb'),'-end');
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'),'-end');
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'),'-end');
cd (laptop_manualROI_folder);
load(template_file);
syllables = [[templates.wavs.segType] -1 100 101]; %[[templates.wavs.segType] 305 -1 100 101]; %
load(annotation_file); 
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
unique_dates = unique(datenum(dates));
cd ROIdata;
%%
% '2017_05_11' '2017_05_15' '2017_05_25' '2017_05_30' ...
%     '2017_06_05' '2017_06_06' '2017_06_08' '2017_06_14' '2017_06_15' ...
%     '2017_06_16' '2017_06_19' '2017_06_20' '2017_06_21' '2017_06_22' '2017_06_23' ...
%     '2017_06_24' '2017_06_25' '2017_06_26' '2017_07_03' '2017_07_04'
SubsetOfDays = {'2017_05_23' '2017_05_25' '2017_05_26'};
%Day = '2017_06_15';
%SubsetOfDays = mat2cell(datestr(unique_dates,'yyyy_mm_dd'),ones(1,numel(unique_dates)),10);
for dayn = 1:numel(SubsetOfDays)
    Day = SubsetOfDays{dayn};
    n_del_frames = 6;
    filecnts = find(datenum(dates) == datenum(Day));
    for cnt = filecnts(1):filecnts(end) %963:numel(ord)
        cd([laptop_manualROI_folder '/ROIdata']);
        hashnum = ord(cnt); %8110;

        loc = min(find(locs == hashnum));
        phrases = return_phrase_times(elements{loc});

        cd(dates(loc,:)); % 2017_07_06;
        try
            load([raw_data_prefix keys{loc}(1:end-3) 'mat']); % lrb85315_8105_2017_07_06_05_37_56;
        catch em
            if strcmp(dates(cnt,:),'2017_04_19')
                continue;
            end
            tokens = regexp(keys{loc},'_','split');
            tokens{2} = num2str(str2num(tokens{2})); %sprintf('%04d',str2num(tokens{2}));
            fname = char(join(tokens,'_'));
            [SUCCESS,MESSAGE,MESSAGEID] = movefile([raw_data_prefix fname(1:end-3) 'mat'],[raw_data_prefix keys{loc}(1:end-3) 'mat']);
            load([raw_data_prefix keys{loc}(1:end-3) 'mat']);
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
end
%cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs/ROIdata');

%% Folders that contain data
% Folders on laptop:
laptop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs';
laptop_wav_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav';
laptop_gif_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/gif';
laptop_storage_folder = '/Volumes/CanaryData/DATA/lrb853_15/movs/';
laptop_annotated_dir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated';
laptop_annotated_images_dir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images';
DamagedFolder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/too_large_or_damaged/';
laptop_manualROI_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs';

%% Single bird
time_padding = 2; %seconds around phrase
bird_name = 'lrb85315';
cd (laptop_manualROI_folder);
load lrb85315template;
syllables = [[templates.wavs.segType] -1 102 103];
n_syllables = numel(syllables);
freq_min = 300; freq_max = 8000;
colors = distinguishable_colors(n_syllables);
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


%% Single day, selected ROIs
Day = '2017_03_29';
sylnum = 5;
ROIs = [5 6]; %18; %[5 6 7]; %[3 9 13];%

warp = 0;
locktoonset = 1;
mulcnt = 3;
opacity_factor = 0.8;

n_del_frames = 4;

h=figure('Visible','on','Position',[77          91        640         600]);
% for i = 1:numel(ROIs)
%     h=figure('Visible','on','Position',[77          91        2215         420]);
%     hs = [hs; h];
% end

cd([laptop_manualROI_folder '/ROIdata/' Day]);

FILES = dir('ROIdata*.mat');
FILES = {FILES.name};
cnt = 0;

for fnum = 1:numel(FILES)
    fname = FILES{fnum};
    tokens = regexp(fname,'_','split');
    loc = find(locs == str2num(tokens{3}));
    phrases = return_phrase_times(elements{loc});
    if ismember(sylnum,phrases.phraseType)
        load(fname);
        display(fname);
        s = zscore(dff(:,n_del_frames+1:end)')';
        t = vidTimes;   
        phrase_locs = find(phrases.phraseType == sylnum);
        for phrase_loc = 1:numel(phrase_locs)
            phrasenum = phrase_locs(phrase_loc);
            tonset = phrases.phraseFileStartTimes(phrasenum);
            toffset = phrases.phraseFileEndTimes(phrasenum);
            for roi_n = 1:numel(ROIs)
%                 if (phrasenum < numel(phrases.phraseFileStartTimes))
%                     subplot(1,numel(ROIs)+1,roi_n);
%                 else
%                     subplot(1,numel(ROIs)+1,numel(ROIs)+1);
%                 end
                subplot(1,numel(ROIs),roi_n);
                %axes(hs(roi_n));
                signal = smooth(s(ROIs(roi_n),:),3);
                timetag = (t(n_del_frames+1:end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
                %signal = (signal-quantile(signal,0.2))/max(signal);
                for currphrase = 1:numel(phrases.phraseType)
                    plot(timetag(t(n_del_frames+1:end) >= phrases.phraseFileStartTimes(currphrase) & ...
                         t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)), ...
                         signal(t(n_del_frames+1:end) >= phrases.phraseFileStartTimes(currphrase) & ...
                         t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase))+cnt*mulcnt, ...
                         'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
                     hold on;
                end
                plot(timetag(t(n_del_frames+1:end) >= phrases.phraseFileEndTimes(currphrase) & ...
                         t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)+2), ...
                         signal(t(n_del_frames+1:end) >= phrases.phraseFileEndTimes(currphrase) & ...
                         t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)+2)+cnt*mulcnt, ...
                         'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','--');
                %plot(,signal+cnt*0.0,'LineWidth',2,'Color',[0 0 0 0.3]);
                cnt = cnt+1;
                set(gca,'color','none');
                box off;
                if (roi_n == 1)
                    title(['Phrase #' num2str(sylnum) ' locked Ca signals from ' datestr(Day,'yyyy-mm-dd')])
                end
                if (roi_n < numel(ROIs))
                    set(gca,'XTick',[]);
                else
                    set(gca,'XTick',[0 1]*locktoonset+(1-locktoonset)*[-1 0]);
                end
                set(gca,'YTick',[]);
                ylabel(['ROI# ' num2str(ROIs(roi_n))]);
                if (roi_n == numel(ROIs))
                    if warp == 1
                        xlabel('Warped Time');
                    else
                        xlabel('Real Time');
                    end
                end
                
                set(gca,'FontSize',16);
                axis tight;
                xlim([-1 3]-(1-locktoonset));
            end
        end
        
    end
end








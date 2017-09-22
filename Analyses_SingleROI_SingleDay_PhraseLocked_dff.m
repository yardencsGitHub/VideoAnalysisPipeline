%% 
% get the structure 'results' that holds the Ca signal variability tests
% (from running 'locate_roi_phrased_locked_activity_types')
% for example:
%load('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs/ROIdata/tmpresults2.mat');
sylables = [0:9 200:209 300:309 400:409 500];
n_syllables = numel(syllables);
thr = 1.2;
bird_name = 'lrb85315';
bird_folder_name = 'lrb853_15';
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];
laptop_manualROI_analyses_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name '/ManualROIs/PhraseLockedSingleDayROIs'];
clear Rs Ps;
for Day_num = 1: size(unique_dates,1)
    Day = unique_dates(Day_num,:);
    nROIs = size(results(Day_num).results_ratio_mean,1);
    for roi_n = 1:nROIs
        for sylnum = 1:n_syllables
            if results(Day_num).results_ratio_mean(roi_n,sylnum) > thr
                [h,r,p] = PhraseLockedSingleDayManualROIs_function(Day,syllables(sylnum),roi_n,1,2);
                Rs(Day_num,roi_n,sylnum) = r;
                Ps(Day_num,roi_n,sylnum) = p;
                fname = [bird_name ' ' datestr(datenum(Day),'yyyy-mm-dd') ' syl=' num2str(syllables(sylnum)) ...
                    ' roi #' num2str(roi_n) ' onst=' num2str(locktoonset) ', spk=' num2str(spikes)];
                set(gca,'CameraPosition',[-11.5361 -174.7864 25.2058]);
                hgsave(h,fullfile(laptop_manualROI_analyses_folder,[fname '.fig']));
                saveas(h,fullfile(laptop_manualROI_analyses_folder,[fname '.png']));
                hgclose(h);
            end
        end
    end
end
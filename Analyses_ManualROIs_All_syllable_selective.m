% Find all cells with activity peaks that exceed double the median among
% syllables. Peaks are definded as the 90% quantile among all repetitions
%%
bird_name = 'lrb85315';
bird_folder_name = 'lrb853_15';
template_file = 'lrb85315template';
annotation_file = 'lrb85315auto_annotation5_fix';

laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

target_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name '/ManualROIs/PhraseLockedSyllablePeaking'];
%% get the 'results' dataset from running 'locate_roi_phrased_locked_activity_types'
highcounts = cell(71,1);
syltypes = {};
daystrings = [];
for dayn = 1:numel(results)
    daystrings = [daystrings; results(dayn).Date];
    tmp = [];
    tmpsyll = {};
    for cellnum = 1:size(results(dayn).Max,1)
        tmp = [tmp sum(results(dayn).Max(cellnum,:) > nanmedian(results(dayn).Max(cellnum,:))*2)];
        tmpsyll = {tmpsyll{:} syllables(results(dayn).Max(cellnum,:) > nanmedian(results(dayn).Max(cellnum,:))*2)};
    end
    highcounts{dayn} = tmp;
    syltypes{dayn} = tmpsyll;
end
%%
cd(laptop_manualROI_folder);
for dayn = 1:numel(results)
    for cellnum = 1:size(results(dayn).Max,1)
        if ~isempty(syltypes{dayn}{cellnum})
            syls = syltypes{dayn}{cellnum};
            for sylnum = 1:numel(syls)
                currsyl = syls(sylnum);
                h=figure('Visible','off'); ax = axes; [ax,r,p] = SequenceLockedSingleDayManualROIs_function(ax,daystrings(dayn,:),currsyl,[nan nan],cellnum,1,2,0,2);
                filename = [bird_name '_' daystrings(dayn,:) '_syl' num2str(currsyl) '_roi_' num2str(cellnum)];
                set(gca,'CameraPosition',[0.4057 -178.3372    4.5232]);
                saveas(h,fullfile(target_dir,[filename '.png']));
                hgsave(h,fullfile(target_dir,filename));
                hgclose(h);
            end
        end
    end
end

%%




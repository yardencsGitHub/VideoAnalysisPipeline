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

target_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name '/ManualROIs/PhraseLockedFWHM'];
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
for dayn = 39:numel(results)
    for cellnum = 1:size(results(dayn).Max,1)
        if ismember(405,syltypes{dayn}{cellnum})
            syltypes{dayn}{cellnum} = unique([setdiff(syltypes{dayn}{cellnum},405) 404]);
        end
        if ~isempty(syltypes{dayn}{cellnum})
            syls = syltypes{dayn}{cellnum};
            for sylnum = 1:numel(syls)
                currsyl = syls(sylnum);
                syllable = currsyl;
                Day = dayn;
                h = figure('Position',[588         234        1898        1104],'Visible','off');
                subplot(3,4,1); plot(results(34).Max(27,:),'ro'); set(gca,'XTick',1:41); set(gca,'XTickLabel',syllables);
                xlabel('syllables'); ylabel('0.9 quantile df/f');
                title(num2str(syltypes{Day}{cellnum}));
                ax = subplot(3,4,2); [ax,r,p] = SequenceLockedSingleDayManualROIs_function(ax,results(Day).Date,syllable,[nan nan],cellnum,1,2,0,2);
                set(gca,'CameraPosition',[-14.5203 -192.6996    0.9790]);
                [xidx,mn,se,mn2,se2,I1,durations] = MeanZscoredDFF_function(results(Day).Date,syllable,cellnum,1,0,0);
                subplot(3,4,6); plot(xidx,nanmean(I1));
                xlabel('Real Time (Sec)'); ylabel('<df/f>')
                activity_max = results(Day).Max(cellnum,syllables==syllable);
                activity_min = results(Day).Min(cellnum,syllables==syllable);
                edges = [0.5 0.5];
                fwhm = [];
                cnts = [];
                for cnt = 1:numel(durations)
                    indx = find(xidx >=0-edges(1) & xidx <= durations(cnt)+edges(2));
                    tmp = intersect(find(I1(cnt,:) > max(I1(cnt,xidx >=0 & xidx <= durations(cnt)))/2),indx);
                    display([cnt xidx([min(tmp) max(tmp)])]);
                    if isempty(tmp)
                        fwhm = [fwhm; nan nan];
                    else
                        fwhm = [fwhm; xidx([min(tmp) max(tmp)])];
                    end
                end

                maxs = [];
                for cnt = 1:numel(durations)
                    indx = find(xidx >=0 & xidx <= durations(cnt));

                    maxs = [maxs; max(I1(cnt,indx))];

                end
                subplot(3,4,5); plot(fwhm(:,:)); hold on; plot(durations(:));
                xlabel('repetition #'); ylabel('Time from onset (sec)'); legend({'early FWHM' 'late FWHM' 'durations'});

                repidx = find(fwhm(:,1)>-0.25 & maxs >= activity_max/2 & durations >= 0.25);
                p = nan, r = p;
                try
                    [r,p] = corr([fwhm(repidx,:) diff(fwhm(repidx,:)')'],durations(repidx)); 
                catch em
                end
                subplot(3,4,9); plot(maxs,'o','MarkerSize',14,'MarkerFaceColor','b'); 
                hold on; plot(setdiff(1:numel(maxs),repidx),maxs(setdiff(1:numel(maxs),repidx)),'o','MarkerSize',14,'MarkerFaceColor',[0.5 0.5 0.5]); 
                line([1 numel(maxs)],[activity_max/2.5 activity_max/2.5],'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','--')
                xlabel('repetition #'); ylabel('df/f'); title(['Maxs ' num2str(r') ', ' num2str(p')]);

                subplot(3,4,10); hist(maxs,20);
                xlabel('Max value'); ylabel('# of repetitions');

                edges = [0.1 0.1];
                fwhm = [];
                cnts = [];
                for cnt = 1:numel(durations)
                    indx = find(xidx >=0-edges(1) & xidx <= durations(cnt)+edges(2));
                    tmp = intersect(find(I1(cnt,:) > max(I1(cnt,xidx >=0 & xidx <= durations(cnt)))/2),indx);
                    display([cnt xidx([min(tmp) max(tmp)])]);
                    if isempty(tmp)
                        fwhm = [fwhm; nan nan];
                    else
                        fwhm = [fwhm; xidx([min(tmp) max(tmp)])];
                    end
                end
                subplot(3,4,7); plot(fwhm(:,:)); hold on; plot(durations(:));
                xlabel('repetition #'); ylabel('Time from onset (sec)'); legend({'early FWHM' 'late FWHM' 'durations'});

                repidx = find(fwhm(:,1)>-0.25 & maxs >= activity_max/2 & durations >= 0.25);
                p = nan, r = p;
                try
                    [r,p] = corr([fwhm(repidx,:) diff(fwhm(repidx,:)')'],durations(repidx)); 
                catch em
                end
                subplot(3,4,11); plot(maxs,'o','MarkerSize',14,'MarkerFaceColor','b'); 
                hold on; plot(setdiff(1:numel(maxs),repidx),maxs(setdiff(1:numel(maxs),repidx)),'o','MarkerSize',14,'MarkerFaceColor',[0.5 0.5 0.5]); 
                line([1 numel(maxs)],[activity_max/2.5 activity_max/2.5],'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','--')
                xlabel('repetition #'); ylabel('df/f'); title(['Maxs ' num2str(r') ', ' num2str(p')]);
                filename = [bird_name '_' daystrings(dayn,:) '_syl' num2str(currsyl) '_roi_' num2str(cellnum)];
                saveas(h,fullfile(target_dir,[filename '.png']));
                hgsave(h,fullfile(target_dir,filename));
                hgclose(h);
            end
        end
    end
end
%%
% figure; plot(xidx,nanmean(I1));
% figure; plot(fwhm(:,:)); hold on; plot(durations(:));

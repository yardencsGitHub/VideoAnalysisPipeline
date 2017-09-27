%%
specificity_thr = 0.5;
%%
bird_name = 'lrb85315';
bird_folder_name = 'lrb853_15';
template_file = 'lrb85315template';
annotation_file = 'lrb85315auto_annotation5';
CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

laptop_manualROI_analyses_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name '/ManualROIs/PhraseLockedSpikeTiming'];
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
%%
syllables = [0:9 200:209 300:309 400:409 500];
long_syllables = [3 6 7 201 205 206 207 307 401 402 403 404 405 406 407 409];
midrange_syllables = [2 200 202 204 300 301 302 303 400 500];
short_syllables = [0 1 4 5 8 9 203 208 209 304 305 306 308 309 408];

cd (laptop_manualROI_folder);

n_syllables = numel(syllables);


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
unique_dates = datestr(setdiff(unique(datenum(dates)),736804),'yyyy_mm_dd'); %does not include 04/19th (remove for other birds)

%% check which rois 'prefer' one or few phrases


%%
if false
    Day_num = 65;
    load(fullfile(laptop_manualROI_analyses_folder,['AlignedSpikes_' bird_name '_' unique_dates(Day_num,:) '.mat']));
    phrase_specificity = [];
    for roi_n = 1:numel(all_syl_num)
        spikes_in_phraseType = cellfun(@numel,{per_syl_num(roi_n,:).data})./(num_phraseType+1e-6)';
        [spikes_in_phraseType, indx] = sort(spikes_in_phraseType,'descend');
        cumm_rate_in_phraseType = cumsum(spikes_in_phraseType)/sum(spikes_in_phraseType);
        phrase_specificity = [phrase_specificity; min(find(cumm_rate_in_phraseType >= specificity_thr))];
        if (phrase_specificity(end) == 1)
            spikes_in_phraseType = cellfun(@numel,{per_syl_num(roi_n,:).data})./(num_phraseType+1e-6)';
            sylnum = syllables(min(find(spikes_in_phraseType == max(spikes_in_phraseType))));
            [h,r,p] = PhraseLockedSingleDayManualROIs_function(unique_dates(Day_num,:),sylnum,roi_n,1,2); set(h,'Visible','on');
            figure; h1 = subplot(1,3,1);
            axes_to_be_copied = findobj(h,'type','axes'); 
            chilred_to_be_copied = get(axes_to_be_copied,'children'); 
            copyobj(chilred_to_be_copied,h1);
            xlim([-1 3]);
            set(gca,'CameraPosition',[-11.5361 -174.7864 25.2058]);
            set(gca,'FontSize',16);
            title(char(h.Children.Title.String));
            hgclose(h);
            subplot(1,3,2); polarhistogram(angle(per_syl_phases(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data),'BinWidth',pi/18);
            set(gca,'FontSize',16);
            subplot(1,3,3); 
            [n,x] = hist(per_syl_time(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data(:,1),0:0.2:2);
            plot(x,n/sum(n),'LineWidth',2,'Color','b'); hold on;
            [n,x] = hist(per_syl_time(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data(:,2),0:0.2:2);
            plot(x,n/sum(n),'LineWidth',2,'Color','r');
            set(gca,'FontSize',16); xlabel('Time (sec) from phrase onset \ offset'); ylabel('Frac. bursts');
            legend({'onset' 'offset'});
            set(gcf,'Position',[352         688        1984         505]);

        end
    end
end
% This can be used to perform histograms

%% repeat. check which roi prefers which syllable but normalize with phrase ratio
% (so it is basically the rate ..)
%sig_type = 1;
zscoring_type = 1;
for Day_num = 21:71
    cd(laptop_manualROI_analyses_folder);
    if ~exist(unique_dates(Day_num,:),'dir')
        mkdir(unique_dates(Day_num,:));
    end
    load(fullfile(laptop_manualROI_analyses_folder,['AlignedSpikes_' bird_name '_' unique_dates(Day_num,:) '.mat']));
    phrase_specificity = [];
    for roi_n = 1:numel(all_syl_num)
        spikes_in_phraseType = cellfun(@numel,{per_syl_num(roi_n,:).data})./(num_phraseType+1e-6)';
        [spikes_in_phraseType, indx] = sort(spikes_in_phraseType,'descend');
        cumm_rate_in_phraseType = cumsum(spikes_in_phraseType)/sum(spikes_in_phraseType);
        phrase_specificity = [phrase_specificity; min(find(cumm_rate_in_phraseType >= specificity_thr))];
        try 
            if (phrase_specificity(end) <= 6  && ~isnan(cumm_rate_in_phraseType(1)))
                for FR_ord = 1:phrase_specificity(end)            
                    %spikes_in_phraseType = cellfun(@numel,{per_syl_num(roi_n,:).data})./(num_phraseType+1e-6)';
                    sylnum = syllables(indx(FR_ord)); %syllables(min(find(spikes_in_phraseType == max(spikes_in_phraseType))));
                    hf = figure('Visible','off'); 

                    h1 = subplot(3,4,[1 2]);
                    [h,r,p] = PhraseLockedSingleDayManualROIs_function(h1,unique_dates(Day_num,:),sylnum,roi_n,1,2,0); %set(h,'Visible','on');      
                    xlim([-1 3]);
                    set(gca,'CameraPosition',[-11.5361 -174.7864 25.2058]);
                    set(gca,'FontSize',16);

                    h2 = subplot(3,4,[3 4]);
                    [h,r,p] = PhraseLockedSingleDayManualROIs_function(h2,unique_dates(Day_num,:),sylnum,roi_n,0,0,0); %set(h,'Visible','on');      
                    xlim([-2 2]);
                    set(gca,'CameraPosition',[-11.5361 -174.7864 25.2058]);
                    set(gca,'FontSize',16);
                    title(['r,p=' num2str(r) ',' num2str(p)]);

                    h3 = subplot(3,4,[5 6]);
                    [h,r,p] = PhraseLockedSingleDayManualROIs_function(h3,unique_dates(Day_num,:),sylnum,roi_n,1,1,1); %set(h,'Visible','on');      
                    xlim([-1 2]);
                    set(gca,'CameraPosition',[-11.5361 -174.7864 25.2058]);
                    set(gca,'FontSize',16);
                    title(['r,p=' num2str(r) ',' num2str(p)]);

                    subplot(3,4,7); polarhistogram(angle(per_syl_phases(roi_n,indx(FR_ord)).data),'BinWidth',pi/18);
                    set(gca,'ThetaTickLabel',[]);
                    set(gca,'FontSize',16);
                    subplot(3,4,8);
                    display([sylnum roi_n])
                    [mn,se] = MeanZscoredDFF_angular_function(unique_dates(Day_num,:),sylnum,roi_n,zscoring_type);
                    polarplot(0:2*pi/1000:2*pi-2*pi/1000,mn,'b','LineWidth',2);
                    hold on
                    polarplot(0:2*pi/1000:2*pi-2*pi/1000,mn+se,'b','LineWidth',1,'LineStyle','--');
                    polarplot(0:2*pi/1000:2*pi-2*pi/1000,mn-se,'b','LineWidth',1,'LineStyle','--');
                    set(gca,'ThetaTickLabel',[]);
                    set(gca,'FontSize',16);
                    subplot(3,4,[9 10]); 
                    [n,x] = hist(per_syl_time(roi_n,indx(FR_ord)).data(:,1),0:0.2:2); %min(find(spikes_in_phraseType == max(spikes_in_phraseType)))
                    plot(x,n/sum(n),'LineWidth',2,'Color','b'); hold on;
                    [n,x] = hist(per_syl_time(roi_n,indx(FR_ord)).data(:,2),0:0.2:2);
                    plot(x,n/sum(n),'LineWidth',2,'Color','r');
                    set(gca,'FontSize',16); xlabel('Time (sec) from phrase edge'); ylabel('Frac. bursts');
                    legend({'onset' 'offset'});
                    title([num2str(numel(per_syl_num(roi_n,indx(FR_ord)).data)) ' spikes, FR=' num2str(spikes_in_phraseType(FR_ord))]);
            %         title([vartest2(per_syl_time(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data(:,1),per_syl_time(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data(:,2),'Tail','right') ...
            %             vartest2(per_syl_time(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data(:,1),per_syl_time(roi_n,min(find(spikes_in_phraseType == max(spikes_in_phraseType)))).data(:,2),'Tail','left')]);
                    [xidx,mn,se,mn2,se2] = MeanZscoredDFF_function(unique_dates(Day_num,:),sylnum,roi_n,1,0,zscoring_type);
                    temp_lcs = ~isnan(mn);
                    subplot(3,4,[11  12]);
                    %fill([xidx(temp_lcs) fliplr(xidx(temp_lcs))],[mn(temp_lcs)+se(temp_lcs) fliplr(mn(temp_lcs)-se(temp_lcs))],'b','FaceAlpha',0.7,'LineStyle','none');
                    plot(xidx(temp_lcs),mn(temp_lcs)+se(temp_lcs),'Color',[0 0 1 0.5],'LineWidth',1);
                    hold on;
                    plot(xidx(temp_lcs),mn(temp_lcs)-se(temp_lcs),'Color',[0 0 1 0.5],'LineWidth',1);
                    plot(xidx(temp_lcs),mn(temp_lcs),'b','LineWidth',2);

                    temp_lcs = ~isnan(mn2);
                    %fill([xidx(temp_lcs) fliplr(xidx(temp_lcs))],[mn2(temp_lcs)+se2(temp_lcs) fliplr(mn2(temp_lcs)-se2(temp_lcs))],'b','FaceAlpha',0.3,'LineStyle','none');
                    plot(xidx(temp_lcs),mn2(temp_lcs)+se2(temp_lcs),'Color',[0 0 1 0.5],'LineWidth',1,'LineStyle','--');
                    hold on;
                    plot(xidx(temp_lcs),mn2(temp_lcs)-se2(temp_lcs),'Color',[0 0 1 0.5],'LineWidth',1,'LineStyle','--');
                    plot(xidx(temp_lcs),mn2(temp_lcs),'b','LineWidth',1,'LineStyle','--');

                    [xidx,mn,se,mn2,se2] = MeanZscoredDFF_function(unique_dates(Day_num,:),sylnum,roi_n,0,0,zscoring_type);
                    temp_lcs = ~isnan(mn);
                    plot(xidx(temp_lcs),mn(temp_lcs)+se(temp_lcs),'Color',[1 0 0 0.5],'LineWidth',1);
                    hold on;
                    plot(xidx(temp_lcs),mn(temp_lcs)-se(temp_lcs),'Color',[1 0 0 0.5],'LineWidth',1);
                    %fill([xidx(temp_lcs) fliplr(xidx(temp_lcs))],[mn(temp_lcs)+se(temp_lcs) fliplr(mn(temp_lcs)-se(temp_lcs))],'r','FaceAlpha',0.7,'LineStyle','none');
                    hold on;
                    plot(xidx(temp_lcs),mn(temp_lcs),'r','LineWidth',2);

                    temp_lcs = ~isnan(mn2);

                    %fill([xidx(temp_lcs) fliplr(xidx(temp_lcs))],[mn2(temp_lcs)+se2(temp_lcs) fliplr(mn2(temp_lcs)-se2(temp_lcs))],'r','FaceAlpha',0.3,'LineStyle','none');
                    plot(xidx(temp_lcs),mn2(temp_lcs)+se2(temp_lcs),'Color',[1 0 0 0.5],'LineWidth',1,'LineStyle','--');
                    hold on;
                    plot(xidx(temp_lcs),mn2(temp_lcs)-se2(temp_lcs),'Color',[1 0 0 0.5],'LineWidth',1,'LineStyle','--');
                    plot(xidx(temp_lcs),mn2(temp_lcs),'r','LineWidth',1,'LineStyle','--');

                    set(gca,'FontSize',16); xlabel('Time (sec) from phrase edge'); ylabel('Mean dff (norm.)');
                    title([FR_ord phrase_specificity(end)])
                    set(gcf,'Position',[98         129        1148        1088]);
                    filename = [bird_name '_' unique_dates(Day_num,:) '_roi' num2str(roi_n) '_syl' num2str(sylnum) ...
                         '_' num2str(phrase_specificity(end)) '_' num2str(FR_ord)];
                    cd(laptop_manualROI_analyses_folder);

                    saveas(hf,fullfile(laptop_manualROI_analyses_folder,unique_dates(Day_num,:),[filename '.png']));
                    hgsave(hf,fullfile(laptop_manualROI_analyses_folder,unique_dates(Day_num,:),[filename '.fig']));
                    hgclose(hf);
                end
            end
        catch eeem
        end
    end
    cd(laptop_manualROI_analyses_folder);
end
%% check transition related activity

function [ax,r,p] = PhraseLockedSingleDayManualROIs_function(ax,Day,sylnum,ROIs,locktoonset,spikes,warp)
%[h,r,p] = PhraseLockedSingleDayManualROIs_function('2017_06_28',301,51,1,2);
%%
thr = 0.02;
zscoring_type = 0;
delete_frames = 1;
n_del_frames = 6;
 
hvc_offset = 0.04;
    %Day = '2017_06_28';
    %sylnum = 200;
    %ROIs = 34; %[304 -->23]; [200 --> 12] %18; %[5 6 7]; %[3 9 13];%
    g = 0.9;
    %warp = 0;
    %locktoonset = 0;
    mulcnt = 2;
    %spikes = 2;
    edges = [0 0];


    opacity_factor = 0.5;
    
   
    %%
    addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
    addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
    bird_name = 'lrb85315';
    bird_folder_name = 'lrb853_15';
    template_file = 'lrb85315template';
    annotation_file = 'lrb85315auto_annotation5';
    CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
    %% Folders that contain data
    % Folders on laptop:
    laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
    laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
    laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
    laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
    laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
    laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
    DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
    laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];
    %% Folders on desktop
    addpath(genpath('/Users/yardenc/Documents/GitHub/small-utils'));
    addpath(genpath(CNMFEfolder));
    %laptop_manualROI_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs';
    %% Single bird
    time_padding = 2; %seconds around phrase
    cd (laptop_manualROI_folder);
    load(template_file);
    syllables = [[templates.wavs.segType] -1 102 103];
    n_syllables = numel(syllables);
    freq_min = 300; freq_max = 8000;
    colors = distinguishable_colors(n_syllables,'w');
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


    %% Single day, selected ROIs

    if isempty(ax)
        h=figure('Visible','off','Position',[77          91        640         600]);
        ax = axes;
    end
    % for i = 1:numel(ROIs)
    %     h=figure('Visible','on','Position',[77          91        2215         420]);
    %     hs = [hs; h];
    % end

    cd([laptop_manualROI_folder '/ROIdata/' Day]);

    FILES = dir('ROIdata*.mat');
    FILES = {FILES.name};
    hits = [];
    durations = [];
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(elements{loc});
        if ismember(sylnum,phrases.phraseType)
            load(fname);
            %s = [s; zscore(dff(ROIs,n_del_frames+1:end)')];

            phrase_locs = find(phrases.phraseType == sylnum);
            for phrase_loc = 1:numel(phrase_locs)
                phrasenum = phrase_locs(phrase_loc);

                tonset = phrases.phraseFileStartTimes(phrasenum);
                toffset = phrases.phraseFileEndTimes(phrasenum);
                hits = [hits; fnum phrasenum];
                durations = [durations; toffset-tonset];
            end
        end
    end
    [durations,dur_idx] = sort(durations);
    hits = hits(dur_idx,:);



    cnt = 0;
    sig_integrals = [];
    sig_std_in = [];
    sig_std_out = [];
    for cnt = 1:size(hits,1)
        fnum = hits(cnt,1);
        phrasenum = hits(cnt,2);
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(elements{loc});
        %if ismember(sylnum,phrases.phraseType)
        load(fname);
        display([num2str(cnt) ' ' fname]);
        
        dff_tmp = dff(:,n_del_frames+1:end);
        y  = dff_tmp;% - ones(size(dff_tmp,1),1)*smooth(mean(dff_tmp),100)';
        if delete_frames == 1
            if zscoring_type == 1
                y = reshape(zscore(y(:)),size(y));                
            end
            t = vidTimes(n_del_frames+1:end)+hvc_offset;
        else
            if zscoring_type == 1
                %
                %y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) y];
                y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) reshape(zscore(y(:)),size(y))];          
                %s = reshape(zscore(y(:)),size(y));
            else
                y = dff;
            end
            t = vidTimes+hvc_offset;
        end
    
        if spikes ==2
            s = y; %[dff(:,1:n_del_frames) y];
            %s = reshape(zscore(y(:)),size(y));
            %s = zscore(detrend(y'))';
        end
        %phrase_locs = find(phrases.phraseType == sylnum);
        %for phrase_loc = 1:numel(phrase_locs)
        %    phrasenum = phrase_locs(phrase_loc);
        tonset = phrases.phraseFileStartTimes(phrasenum);
        toffset = phrases.phraseFileEndTimes(phrasenum);
        for roi_n = 1:numel(ROIs)
            if spikes < 2
                try %detrend
                    %[c, s, options] = deconvolveCa(detrend(y(roi_n,:)),'ar1',g,'method','constrained-foopsi');
                    %[c, s, options] =
                    %deconvolveCa((y(ROIs(roi_n),:)),'ar2','method','thresholded','optimize_b',1); [1.3 -0.422]
                    [c, s, options] = deconvolveCa(y(ROIs(roi_n),:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');%,'optimize_pars');
                catch em
                    [c, s, options] = deconvolveCa((y(ROIs(roi_n),:)),'ar2','method','foopsi','optimize_b',1);
                    %[c, s, options] = deconvolveCa(detrend(y(roi_n,:)),'ar1',g,'method','foopsi');
                end
                %[c, s, options] = deconvolveCa((y(ROIs(roi_n),n_del_frames+1:end)),'ar1',g,'method','constrained-foopsi');  
            end
            switch spikes
                case 1
                    signal = s;
                case 0
                  signal = c; %smooth(s(ROIs(roi_n),:),3);  
                otherwise
                    signal = smooth(s(ROIs(roi_n),:),3);
            end
            %n_del_frames = 0;
            sig_integrals = [sig_integrals; ...
                sum(signal((t >= tonset) & (t <= toffset)))];
            sig_std_in = [sig_std_in; ...
                quantile(signal((t >= (tonset-edges(1))) & (t <= (toffset+edges(2)))),0.9) - ...
                quantile(signal((t >= (tonset-edges(1))) & (t <= (toffset+edges(2)))),0.1)];
            sig_std_out = [sig_std_out; ...
                quantile(signal((t <= (tonset-edges(1))) | (t >= (toffset+edges(2)))),0.9) - ...
                quantile(signal((t <= (tonset-edges(1))) | (t >= (toffset+edges(2)))),0.1)];
            timetag = (t-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
            %signal = (signal-quantile(signal,0.2))/max(signal);
            for currphrase = 1:numel(phrases.phraseType)
                plot3(ax,timetag(t >= phrases.phraseFileStartTimes(currphrase) & ...
                     t <= phrases.phraseFileEndTimes(currphrase)), ...
                     cnt*mulcnt*ones(1,sum(t >= phrases.phraseFileStartTimes(currphrase) & ...
                     t <= phrases.phraseFileEndTimes(currphrase))), ...
                     signal(t >= phrases.phraseFileStartTimes(currphrase) & ...
                     t <= phrases.phraseFileEndTimes(currphrase))+0, ...
                     'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
                 hold on;
                 if (currphrase < numel(phrases.phraseType))
                     startidx = max(find(t <= phrases.phraseFileEndTimes(currphrase)));
                     stopidx = min(find(t >= phrases.phraseFileStartTimes(currphrase+1)));
                     plot3(ax,timetag(startidx:stopidx), ...
                     cnt*mulcnt*ones(1,numel(startidx:stopidx)), ...
                     signal(startidx:stopidx)+0, ...
                     'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','-');
                 end
            end
            plot3(ax,timetag(t >= phrases.phraseFileEndTimes(currphrase) & ...
                     t <= phrases.phraseFileEndTimes(currphrase)+2), ...
                     cnt*mulcnt*ones(1,sum(t >= phrases.phraseFileEndTimes(currphrase) & ...
                     t <= phrases.phraseFileEndTimes(currphrase)+2)), ...
                     signal(t >= phrases.phraseFileEndTimes(currphrase) & ...
                     t <= phrases.phraseFileEndTimes(currphrase)+2)+0, ...
                     'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','--');
            %plot(,signal+cnt*0.0,'LineWidth',2,'Color',[0 0 0 0.3]);
            %cnt = cnt+1;
            set(gca,'color','none');
            %box off;
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
        %end

        %end
    end
    [r, p] = corr(durations,sig_integrals);
    title([bird_name ', ' datestr(datenum(Day),'yyyy-mm-dd') ', syl=' num2str(sylnum) ', roi #' num2str(ROIs) ', corr(r,p)=(' sprintf('%2.2f',r) ...
        ',' sprintf('%1.4f',p) '), onst=' num2str(locktoonset) ', spk=' num2str(spikes)]);
end






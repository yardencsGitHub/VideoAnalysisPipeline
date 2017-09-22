function [xidx,mn,se,mn2,se2] = MeanZscoredDFF_function(Day,sylnum,ROIs,locktoonset,warp,zscoring_type)
%%
    %zscoring_type = 0;
    delete_frames = 1;
    n_del_frames = 6;
    hvc_offset = 0.035;
%% Folders that contain data
    % Folders on laptop:
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
    %% Folders on desktop:
    %addpath(genpath('/Users/yardenc/Documents/GitHub/small-utils'));

    %laptop_manualROI_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs';
    %% Single bird
    time_padding = 2; %seconds around phrase

    cd (laptop_manualROI_folder);
    load lrb85315template;
    syllables = [[templates.wavs.segType] -1 102 103];
    n_syllables = numel(syllables);
    freq_min = 300; freq_max = 8000;
    colors = distinguishable_colors(n_syllables,'k');
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
    %Day = '2017_06_28';
    %sylnum = 5;
    %ROIs = [44]; %18; %[5 6 7]; %[3 9 13];%
   
    %warp = 0;
    %locktoonset = 1;
    %bg_lim = 0.5;

    cd([laptop_manualROI_folder '/ROIdata/' Day]);

    FILES = dir('ROIdata*.mat');
    FILES = {FILES.name};
    %%
    hits = [];
    durations = [];
    s = [];
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(elements{loc});
        if ismember(sylnum,phrases.phraseType)
            load(fname);
            s = [s; zscore(dff(ROIs,n_del_frames+1:end)')];

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
    maxs = max(s);
    [durations,dur_idx] = sort(durations);
    hits = hits(dur_idx,:);
    %%

    if warp == 0
        I1 = zeros(size(hits,1),ceil(max(durations)*1000+4200)); 
    else
        I1 = zeros(size(hits,1),ceil(max(durations)/min(durations)*1000+4200));
    end
    I1 = I1*nan;
    I2 = I1;
    for cnt = 1:size(hits,1)
        fnum = hits(cnt,1);
        phrasenum = hits(cnt,2);
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(elements{loc});
        load(fname);
        display(fname);
        
        dff_tmp = dff(:,n_del_frames+1:end);
        y  = dff_tmp;% - ones(size(dff_tmp,1),1)*smooth(mean(dff_tmp),100)';
        if delete_frames == 1
            if zscoring_type == 1
                s = reshape(zscore(y(:)),size(y));                
            else
                s = y;
            end
            t = vidTimes(n_del_frames+1:end)-hvc_offset;
        else
            if zscoring_type == 1
                %
                %y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) y];
                s = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) reshape(zscore(y(:)),size(y))];          
                %s = reshape(zscore(y(:)),size(y));
            else
                s = dff;
            end
            t = vidTimes-hvc_offset;
        end
        %dff = dff(:,n_del_frames+1:end);
        %y  = dff - ones(size(dff,1),1)*smooth(mean(dff),100)';
        
        
        %
        %t = vidTimes;  
        %[c, s1, options] = deconvolveCa(detrend(s(ROIs,:)),'ar1',0.99,'method','constrained-foopsi');
        tonset = phrases.phraseFileStartTimes(phrasenum);
        toffset = phrases.phraseFileEndTimes(phrasenum);
        signal = (s(ROIs,:)); %smooth(s(ROIs,:),3);
        %signal(signal < quantile(signal,bg_lim)) = 0;
        %signal = signal / max(signal); %maxs
        %n_del_frames = 0;
        %timetag = (t(n_del_frames+1:end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        timetag = (t-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        
        durfactor = 1;
        %if warp == 0
            signal = interp1(timetag,signal,timetag(1):0.001:timetag(end)+1/30);
            %t_temp = interp1(timetag,t_temp,timetag(1):0.001:timetag(end)+1/30);
            timetag = interp1(timetag,timetag,timetag(1):0.001:timetag(end)+1/30);

        idxmap = [1:numel(timetag)] - min(find(abs(timetag) == min(abs(timetag))))+2100+round((1-locktoonset)*max(durations)*(1000*(1-warp)+warp*durfactor));


        t_on = (phrases.phraseFileStartTimes(1)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        t_off = (phrases.phraseFileEndTimes(end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        idxs = find((timetag >= t_on) & (timetag <= t_off) & ...
                (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
                (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
        I1(cnt,idxmap(idxs)) = signal(idxs);

        t_on = (phrases.phraseFileEndTimes(end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        t_off = (phrases.phraseFileEndTimes(end)+2-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        idxs = find((timetag >= t_on) & (timetag <= t_off) & ...
                (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
                (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
        if ~isempty(idxs)
            I1(cnt,idxmap(idxs)) = signal(idxs);     
        end
        I2(cnt,:) = I1(cnt,:);
        t_on = (phrases.phraseFileStartTimes(1)-2-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        t_off = (phrases.phraseFileEndTimes(end)+2-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        t_onset = (tonset-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        t_offset = (toffset-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        idxs = find((timetag < t_onset) & (timetag > t_on) & ...
                (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
                (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
        I2(cnt,idxmap(idxs)) = nan;
        idxs = find((timetag > t_offset) & (timetag < t_off) & ...
                (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
                (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
        I2(cnt,idxmap(idxs)) = nan;
    end

    %%
    %h2 = figure;
    mn = nanmean(I1);
    mn2 = nanmean(I2);
    se2 = nanstd(I2)/sqrt(size(I1,1));
    se = nanstd(I1)/sqrt(size(I1,1));
    xidx = -2.1-max(durations)*(1-locktoonset):1/1000:2.1+max(durations)*(locktoonset);
    if (numel(xidx) > numel(mn))
        xidx = xidx(1:numel(mn));
    end
     if (numel(xidx) < numel(mn))
         mn = mn(1:numel(xidx));
         se = se(1:numel(xidx));
     end
    % fill([xidx fliplr(xidx)],[mn+se fliplr(mn-se)],[1 0 0],'FaceAlpha',0.5);
    % hold on;
    % plot(-2.1-max(durations)*(1-locktoonset):1/1000:2.1+max(durations)*(locktoonset),mean(I1),'b','LineWidth',2);
end








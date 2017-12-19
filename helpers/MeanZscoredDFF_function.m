function [xidx,mn,se,mn2,se2,I1,durations] = MeanZscoredDFF_function(Day,sylnum,flanking_syllabels,ROIs,locktoonset,warp,zscoring_type,spikes)
%%
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation4'};
bird_params = bird1_params;

%     bird_name = 'lbr3022'; %'lrb85315';
%     bird_folder_name = 'lbr3022'; % 'lrb853_15';
%     template_file = 'lbr3022_template.mat'; %'lrb85315template';
%     annotation_file = 'lrb85315auto_annotation5_fix'; %'lbr3022auto_annotation4.mat'; %
    
    bird_name = bird_params{1}; % 'lrb85315'; %'lbr3022'; %
    bird_folder_name = bird_params{2}; %'lrb853_15'; %'lbr3022'; %
    template_file = bird_params{3}; %'lrb85315template'; %'lbr3022_template';%
    annotation_file = bird_params{4};
    %zscoring_type = 0;
    delete_frames = 1;
    n_del_frames = 5;
    hvc_offset = 0.04;
    combine407 = 0;
    remove_hits = [];
    syls_to_combine = [405 404];
    plane_flag = 0;
    align_syls = 0;
%% Folders that contain data
    % Folders on laptop:
    
   

% Folders on laptop:

    
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
    load(template_file);% lrb85315template;
    syllables = [[templates.wavs.segType] -1]; % 102 103];
    n_syllables = numel(syllables);
    freq_min = 300; freq_max = 8000;
    colors = distinguishable_colors(n_syllables,'k');
    load(annotation_file);% lrb85315auto_annotation5_fix;
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
    hits = [];
    durations = [];
    dur_flanking = [];
    max_syls = 0;
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        if combine407
            elements{loc}.segType(elements{loc}.segType==syls_to_combine(1))=syls_to_combine(2);
        end
        phrases = return_phrase_times(elements{loc});
        
        if ismember(sylnum,phrases.phraseType)
            load(fname);
            %s = [s; zscore(dff(ROIs,n_del_frames+1:end)')];

            phrase_locs = find(phrases.phraseType == sylnum);
            for phrase_loc = 1:numel(phrase_locs)
                phrasenum = phrase_locs(phrase_loc);

                tonset = phrases.phraseFileStartTimes(phrasenum);
                toffset = phrases.phraseFileEndTimes(phrasenum);
                syls_in_phrase = find(elements{loc}.segFileEndTimes >= tonset & ...
                elements{loc}.segFileStartTimes <= toffset);
                
                if phrasenum == 1
                    phst(1)=-1;
                else
                    phst(1) = phrases.phraseType(phrasenum-1);
                end
                for num_future = 1:(numel(flanking_syllabels)-1)
                    if phrasenum >= numel(phrases.phraseType)-num_future+1
                        phst(num_future+1)=-1;
                    else
                        phst(num_future+1) = phrases.phraseType(phrasenum+num_future);
                    end
                end
                %phrase_strcture = [phrase_strcture; phst];
                if (nansum(abs(phst-flanking_syllabels)) == 0)
                    max_syls = max(max_syls,numel(syls_in_phrase));
                    hits = [hits; fnum phrasenum phst];
                    durations = [durations; toffset-tonset];
                    
                    if phst(1) == -1
                        drf = [0];
                    else
                        drf = [phrases.phraseFileEndTimes(phrasenum-1)-phrases.phraseFileStartTimes(phrasenum-1)];
                    end
                    if phst(2) == -1
                        drf = [drf 0];
                    else
                        drf = [drf phrases.phraseFileEndTimes(phrasenum+1)-phrases.phraseFileStartTimes(phrasenum+1)];
                    end
                    dur_flanking = [dur_flanking; drf];
                end
            end
        end
    end
   
    [durations,dur_idx] = sort(durations);
    hits = hits(dur_idx,:);
    
    if ~isempty(remove_hits)
        hits(remove_hits,:)=[];
        durations(remove_hits)=[];
    end
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
            t = vidTimes(n_del_frames+1:end)+hvc_offset;
        else
            if zscoring_type == 1
                %
                %y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) y];
                s = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) reshape(zscore(y(:)),size(y))];          
                %s = reshape(zscore(y(:)),size(y));
            else
                s = dff;
            end
            t = vidTimes+hvc_offset;
        end
       
        tonset = phrases.phraseFileStartTimes(phrasenum);
        toffset = phrases.phraseFileEndTimes(phrasenum);
        
        %signal = smooth(s(ROIs,:),3); %smooth(s(ROIs,:),3);
        if spikes < 2
            try %detrend
                %[c, s, options] = deconvolveCa(detrend(y(roi_n,:)),'ar1',g,'method','constrained-foopsi');
                %[c, s, options] =
                %deconvolveCa((y(ROIs(roi_n),:)),'ar2','method','thresholded','optimize_b',1); [1.3 -0.422]
                [c, s, options] = deconvolveCa(y(ROIs,:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');%,'optimize_pars');
            catch em
                [c, s, options] = deconvolveCa((y(ROIs,:)),'ar2','method','foopsi','optimize_b',1);
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
                signal = smooth(s(ROIs,:),3);
        end
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








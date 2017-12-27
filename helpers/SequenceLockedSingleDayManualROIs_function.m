function [ax,r,p] = SequenceLockedSingleDayManualROIs_function(ax,Day,sylnum,flanking_syllabels,ROIs,locktoonset,spikes,warp,order_flag)
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation4'};
bird_params = bird1_params;
% Inputs:
% ax - axes where to plot
% Day - text rep. of date
% sylnum - tag of syllable
% flanking_syllables - 2 vector of tags of syllables before and after
% ROIs - the ROI in question
%%
combine407 = 0;
remove_hits = [];
syls_to_combine = [405 404];
plane_flag = 0;
align_syls = 0;
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
    edges = [0.0 0.0];
    opacity_factor = 0.5;
    
   
    %%
    addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
    addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
    bird_name = bird_params{1}; % 'lrb85315'; %'lbr3022'; %
    bird_folder_name = bird_params{2}; %'lrb853_15'; %'lbr3022'; %
    template_file = bird_params{3}; %'lrb85315template'; %'lbr3022_template';%
    annotation_file = bird_params{4}; %'lrb85315auto_annotation5_fix'; %'lbr3022auto_annotation4';% 
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
    load('syl_dur_gap_stat.mat');
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

    FILES = dir('NonoverlapBaseROIdata*.mat');
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
    
    switch order_flag
        case 1
            [durations,dur_idx] = sort(dur_flanking(:,1)); %dur_flanking(:,2)
        case 2
            [durations,dur_idx] = sort(durations);
        case 3
            [durations,dur_idx] = sort(dur_flanking(:,2));
        case 12
            [durations,dur_idx] = sort(dur_flanking(:,1)+durations);
        case 23
            [durations,dur_idx] = sort(dur_flanking(:,2)+durations);
        case 123
            [durations,dur_idx] = sort(dur_flanking(:,1)+dur_flanking(:,2)+durations);
        case -1
            [durations,dur_idx] = sort(hits(:,3));
        case -3
            [durations,dur_idx] = sort(hits(:,4));
        otherwise
            dur_idx = 1:numel(durations);
    end
     hits = hits(dur_idx,:);

    if ~isempty(remove_hits)
        hits(remove_hits,:)=[];
        durations(remove_hits)=[];
    end

    cnt = 0;
    sig_integrals = [];
    sig_std_in = [];
    sig_std_out = [];
    max_signal = 0;
    for cnt = 1:size(hits,1)
        fnum = hits(cnt,1);
        phrasenum = hits(cnt,2);
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        if combine407
            %elements{loc}.segType(elements{loc}.segType==405)=404;
            elements{loc}.segType(elements{loc}.segType==syls_to_combine(1))=syls_to_combine(2);
        end
        phrases = return_phrase_times(elements{loc});
        %if ismember(sylnum,phrases.phraseType)
        load(fname);
        %display([num2str(cnt) ' ' fname]);
        
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
        display([num2str([cnt tonset]) ' ' fname]);
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
                case 3
%                     events_fname = ['Events_' FILES{fnum}];
%                     load(events_fname);
%                     signal = map_states(ROIs(roi_n),n_del_frames+1:end);
                    sig = y(ROIs(roi_n),:);
                    clear sigma;
                    nstates = 2;
                    dq = (0.5-0.05)/(nstates-2+1e-5);          
                    mu = [quantile(sig,0.05+[0:nstates-2]*dq) max(sig) ];
                    sigma(1,1,:) = [0.05*ones(1,nstates)];
                    CPD = condGaussCpdCreate(mu, sigma);
                    [model, loglikHist] = hmmFit(sig, nstates, 'gauss','pi0',[zeros(1,nstates-1) 1],'emission0',CPD);%,'transPrior',[10 1;10 3],'piPrior',[100 1]);                    
                    path = hmmMap(model, sig)-1;
                  signal = abs(median(path)-path);
                 otherwise
                     signal = smooth(s(ROIs(roi_n),:),3);
            end
            %n_del_frames = 0;
           
            sig_integrals = [sig_integrals; ...
                sum(signal((t >= tonset-edges(1)) & (t <= toffset+edges(2))))];
%             sig_integrals = [sig_integrals; ...
%                 max(signal((t >= tonset-edges(1)) & (t <= tonset+edges(2))))];
            sig_std_in = [sig_std_in; ...
                quantile(signal((t >= (tonset-edges(1))) & (t <= (toffset+edges(2)))),0.9) - ...
                quantile(signal((t >= (tonset-edges(1))) & (t <= (toffset+edges(2)))),0.1)];
            sig_std_out = [sig_std_out; ...
                quantile(signal((t <= (tonset-edges(1))) | (t >= (toffset+edges(2)))),0.9) - ...
                quantile(signal((t <= (tonset-edges(1))) | (t >= (toffset+edges(2)))),0.1)];
            timetag = (t-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
            
            if align_syls == 1
                msyl = m_syls(find(syllables == phrases.phraseType(phrasenum)))/((toffset-tonset)*warp+1-warp);
                mgap = m_gaps(find(syllables == phrases.phraseType(phrasenum)))/((toffset-tonset)*warp+1-warp);
                syls_in_phrase = find(elements{loc}.segFileEndTimes >= phrases.phraseFileStartTimes(phrasenum) & ...
                elements{loc}.segFileStartTimes <= phrases.phraseFileEndTimes(phrasenum));
                n_syls_in_phrase = numel(syls_in_phrase);
                syl_start_times = (elements{loc}.segFileStartTimes(syls_in_phrase)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
                syl_end_times = (elements{loc}.segFileEndTimes(syls_in_phrase)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
                X = [syl_start_times syl_end_times];
                V = tonset + [0:(msyl+mgap):(n_syls_in_phrase-1)*(msyl+mgap) ...
                    (0:(msyl+mgap):(n_syls_in_phrase-1)*(msyl+mgap))+msyl];
                V(end) = toffset;
                V = (V-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
                signal = interp1(t,signal,t(1):0.001:t(end)+1/30);
                timetag = interp1(t,timetag,t(1):0.001:t(end)+1/30);
                timetag(isnan(timetag)) = 0;
                t = interp1(t,t,t(1):0.001:t(end)+1/30);
                tq = interp1(X,V,timetag(t >= phrases.phraseFileStartTimes(phrasenum) & ...
                     t <= phrases.phraseFileEndTimes(phrasenum)));
                 timetag(t >= phrases.phraseFileStartTimes(phrasenum) & ...
                     t <= phrases.phraseFileEndTimes(phrasenum)) = tq;
                max_signal = max(max_signal,max(signal(t >= phrases.phraseFileStartTimes(phrasenum) & ...
                     t <= phrases.phraseFileEndTimes(phrasenum))));
                 if (max_syls <= n_syls_in_phrase)
                     Vsave = V;
                 end
                
            end
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
    
    
    if (plane_flag == 1)           
        for slcnt = 1:max_syls
           plot3(ax,ones(1,size(hits,1))*Vsave(slcnt),[1:size(hits,1)]*mulcnt,ones(1,size(hits,1))*0.0,'Color','g','LineStyle','--');
           plot3(ax,ones(1,size(hits,1))*Vsave(slcnt),ones(1,size(hits,1))*size(hits,1)*mulcnt,[1:size(hits,1)]*max_signal/size(hits,1),'Color','g','LineStyle','--');
           plot3(ax,ones(1,size(hits,1))*Vsave(max_syls+slcnt),[1:size(hits,1)]*mulcnt,ones(1,size(hits,1))*0.0,'Color','r','LineStyle','--');
           plot3(ax,ones(1,size(hits,1))*Vsave(max_syls+slcnt),ones(1,size(hits,1))*size(hits,1)*mulcnt,[1:size(hits,1)]*max_signal/size(hits,1),'Color','r','LineStyle','--');
        end
    end
     [r, p] = corr(durations,sig_integrals);     
    if order_flag < 0
        [p,ANOVATAB,STATS] = anova1(sig_integrals,durations);
        
    else
        figure; plot(durations,sig_integrals,'bo','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','none');
        set(gca,'FontSize',16); xlabel('Durations'); ylabel('Signal Integral')
    end
   
    title([bird_name ', ' datestr(datenum(Day),'yyyy-mm-dd') ', syl=' num2str(sylnum) ', roi #' num2str(ROIs) ', corr(r,p)=(' sprintf('%2.2f',r) ...
      ',' sprintf('%1.4f',p) '), onst=' num2str(locktoonset) ', spk=' num2str(spikes)]);
end






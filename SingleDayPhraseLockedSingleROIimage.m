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
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
%% Folders on desktop:
%addpath(genpath('/Users/yardenc/Documents/GitHub/small-utils'));

%laptop_manualROI_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs';
%% Single bird
time_padding = 2; %seconds around phrase
bird_name = 'lrb85315';
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
Day = '2017_06_28';
sylnum = 5;
ROIs = [44]; %18; %[5 6 7]; %[3 9 13];%
zscoring_type = 1;
warp = 0;
locktoonset = 1;
bg_lim = 0.5;

n_del_frames = 6;
hvc_offset = 0.035;
h=figure('Visible','on','Position',[77          91        640         600]);
% for i = 1:numel(ROIs)
%     h=figure('Visible','on','Position',[77          91        2215         420]);
%     hs = [hs; h];
% end

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
    I = zeros(size(hits,1),ceil(max(durations)*1000+4200),3);
    I1 = zeros(size(hits,1),ceil(max(durations)*1000+4200));
else
    I = zeros(size(hits,1),ceil(max(durations)/min(durations)*1000+4200),3);
    I1 = zeros(size(hits,1),ceil(max(durations)/min(durations)*1000+4200));
end
for cnt = 1:size(hits,1)
    fnum = hits(cnt,1);
    phrasenum = hits(cnt,2);
    fname = FILES{fnum};
    tokens = regexp(fname,'_','split');
    loc = find(locs == str2num(tokens{3}));
    phrases = return_phrase_times(elements{loc});
    load(fname);
    display(fname);
    dff = dff(:,n_del_frames+1:end);
    y  = dff - ones(size(dff,1),1)*smooth(mean(dff),100)';
    if zscoring_type == 1
        s = reshape(zscore(y(:)),size(y));
    else
        s = zscore(y')';
    end
    t = vidTimes+hvc_offset;
    %
    %t = vidTimes;  
    tonset = phrases.phraseFileStartTimes(phrasenum);
    toffset = phrases.phraseFileEndTimes(phrasenum);
    signal = smooth(s(ROIs,:),3);
    signal(signal < quantile(signal,bg_lim)) = 0;
    signal = signal / max(signal); %maxs
    timetag = (t(n_del_frames+1:end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
    t_temp = t(n_del_frames+1:end);
    durfactor = 1;
    %if warp == 0
        signal = interp1(timetag,signal,timetag(1):0.001:timetag(end)+1/30);
        t_temp = interp1(timetag,t_temp,timetag(1):0.001:timetag(end)+1/30);
        timetag = interp1(timetag,timetag,timetag(1):0.001:timetag(end)+1/30);
%     else
%         durfactor = max(durations)/durations(cnt);
%         signal = interp1(timetag,signal,timetag(1):1/durfactor/1000:timetag(end)+1/30);
%         t_temp = interp1(timetag,t_temp,timetag(1):1/durfactor/1000:timetag(end)+1/30);
%         timetag = interp1(timetag,timetag,timetag(1):1/durfactor/1000:timetag(end)+1/30);
%     end
    idxmap = [1:numel(timetag)] - min(find(abs(timetag) == min(abs(timetag))))+2100+round((1-locktoonset)*max(durations)*(1000*(1-warp)+warp*durfactor));
    
    
    t_on = (phrases.phraseFileStartTimes(1)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
    t_off = (phrases.phraseFileEndTimes(end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
    idxs = find((timetag >= t_on) & (timetag <= t_off) & ...
            (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
            (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
    I1(cnt,idxmap(idxs)) = signal(idxs);
    for currphrase = 1:numel(phrases.phraseType)
        t_on = (phrases.phraseFileStartTimes(currphrase)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        t_off = (phrases.phraseFileEndTimes(currphrase)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
        idxs = find((timetag >= t_on) & (timetag <= t_off) & ...
            (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
            (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
        if ~isempty(idxs)
           %plot(timetag(idxs),signal(idxs),'Color',colors(find(syllables == phrases.phraseType(currphrase)),:));
           % hold on;
            I(cnt,idxmap(idxs),:) = signal(idxs)'*colors(find(syllables == phrases.phraseType(currphrase)),:);% + ...
                %(1-signal(idxs)')*(1-colors(find(syllables == phrases.phraseType(currphrase)),:));
            
        end
    end
    t_on = (phrases.phraseFileEndTimes(currphrase)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
    t_off = (phrases.phraseFileEndTimes(currphrase)+2-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
    idxs = find((timetag >= t_on) & (timetag <= t_off) & ...
            (timetag >= -2*locktoonset  + (-max(durations)*(1-warp)-2-warp)*(1-locktoonset)) & ...
            (timetag <= 2*(1-locktoonset)+locktoonset*(2+max(durations)*(1-warp)+warp)));
    if ~isempty(idxs)
        %plot(timetag(idxs),signal(idxs),'Color',colors(find(syllables == phrases.phraseType(currphrase)),:));
        %hold on;
        I(cnt,idxmap(idxs),:) = signal(idxs)'*[0.5 0.5 0.5];% + ...
            %(1-signal(idxs)')*(1-colors(find(syllables == phrases.phraseType(currphrase)),:));
        I1(cnt,idxmap(idxs)) = signal(idxs);
    end
end

%%
image([-2.1-max(durations)*(1-locktoonset) 2.1+max(durations)*(locktoonset)],[1 cnt],I);
if warp
    image(([1 size(I,2)]-2100)/1000,[1 cnt],I);
end

title(['Phrase #' num2str(sylnum) ' locked Ca signals from ' datestr(Day,'yyyy-mm-dd')])
set(gca,'XTick',[0 1]*locktoonset+(1-locktoonset)*[-1 0]);
set(gca,'YTick',[]);
ylabel(['ROI# ' num2str(ROIs)]);

    if warp == 1
        xlabel('Warped Time');
    else
        xlabel('Real Time');
    end


set(gca,'FontSize',16);
%%
h = figure;
mn = mean(I1);
se = std(I1)/sqrt(size(I1,1));
xidx = -2.1-max(durations)*(1-locktoonset):1/1000:2.1+max(durations)*(locktoonset);
fill([xidx fliplr(xidx)],[mn+se fliplr(mn-se)],[1 0 0],'FaceAlpha',0.5);
hold on;
plot(-2.1-max(durations)*(1-locktoonset):1/1000:2.1+max(durations)*(locktoonset),mean(I1),'b','LineWidth',2);




%         plot(timetag(t(n_del_frames+1:end) >= phrases.phraseFileStartTimes(currphrase) & ...
%              t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)), ...
%              signal(t(n_del_frames+1:end) >= phrases.phraseFileStartTimes(currphrase) & ...
%              t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase))+cnt*mulcnt, ...
%              'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
%          hold on;
%    
%     plot(timetag(t(n_del_frames+1:end) >= phrases.phraseFileEndTimes(currphrase) & ...
%              t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)+2), ...
%              signal(t(n_del_frames+1:end) >= phrases.phraseFileEndTimes(currphrase) & ...
%              t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)+2)+cnt*mulcnt, ...
%              'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','--');
%     
%     
% 
% %%
% for fnum = 1:numel(FILES)
%     fname = FILES{fnum};
%     tokens = regexp(fname,'_','split');
%     loc = find(locs == str2num(tokens{3}));
%     phrases = return_phrase_times(elements{loc});
%     if ismember(sylnum,phrases.phraseType)
%         load(fname);
%         display(fname);
%         s = zscore(dff(:,n_del_frames+1:end)')';
%         t = vidTimes;   
%         phrase_locs = find(phrases.phraseType == sylnum);
%         for phrase_loc = 1:numel(phrase_locs)
%             phrasenum = phrase_locs(phrase_loc);
%             tonset = phrases.phraseFileStartTimes(phrasenum);
%             toffset = phrases.phraseFileEndTimes(phrasenum);
%             for roi_n = 1:numel(ROIs)
% %                 if (phrasenum < numel(phrases.phraseFileStartTimes))
% %                     subplot(1,numel(ROIs)+1,roi_n);
% %                 else
% %                     subplot(1,numel(ROIs)+1,numel(ROIs)+1);
% %                 end
%                 subplot(1,numel(ROIs),roi_n);
%                 %axes(hs(roi_n));
%                 signal = smooth(s(ROIs(roi_n),:),3);
%                 timetag = (t(n_del_frames+1:end)-tonset*locktoonset-(1-locktoonset)*toffset)/((toffset-tonset)*warp+1-warp);
%                 %signal = (signal-quantile(signal,0.2))/max(signal);
%                 for currphrase = 1:numel(phrases.phraseType)
%                     plot(timetag(t(n_del_frames+1:end) >= phrases.phraseFileStartTimes(currphrase) & ...
%                          t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)), ...
%                          signal(t(n_del_frames+1:end) >= phrases.phraseFileStartTimes(currphrase) & ...
%                          t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase))+cnt*mulcnt, ...
%                          'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
%                      hold on;
%                 end
%                 plot(timetag(t(n_del_frames+1:end) >= phrases.phraseFileEndTimes(currphrase) & ...
%                          t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)+2), ...
%                          signal(t(n_del_frames+1:end) >= phrases.phraseFileEndTimes(currphrase) & ...
%                          t(n_del_frames+1:end) <= phrases.phraseFileEndTimes(currphrase)+2)+cnt*mulcnt, ...
%                          'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','--');
%                 %plot(,signal+cnt*0.0,'LineWidth',2,'Color',[0 0 0 0.3]);
%                 cnt = cnt+1;
%                 set(gca,'color','none');
%                 box off;
%                 if (roi_n == 1)
%                     title(['Phrase #' num2str(sylnum) ' locked Ca signals from ' datestr(Day,'yyyy-mm-dd')])
%                 end
%                 if (roi_n < numel(ROIs))
%                     set(gca,'XTick',[]);
%                 else
%                     set(gca,'XTick',[0 1]*locktoonset+(1-locktoonset)*[-1 0]);
%                 end
%                 set(gca,'YTick',[]);
%                 ylabel(['ROI# ' num2str(ROIs(roi_n))]);
%                 if (roi_n == numel(ROIs))
%                     if warp == 1
%                         xlabel('Warped Time');
%                     else
%                         xlabel('Real Time');
%                     end
%                 end
%                 
%                 set(gca,'FontSize',16);
%                 axis tight;
%                 xlim([-1 3]-(1-locktoonset));
%             end
%         end
%         
%     end
% end








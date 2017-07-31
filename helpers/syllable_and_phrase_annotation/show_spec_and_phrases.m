%%
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/zftftb'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images';
template_file = 'lrb85315template.mat';
cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
load lrb85315auto_annotation4;
load(template_file);
syllables = [[templates.wavs.segType] -1 102 103];
cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/mat');
n_syllables = numel(syllables);
freq_min = 300; freq_max = 8000;
colors = distinguishable_colors(n_syllables);
%%
<<<<<<< Updated upstream
for fnum = 2669:2669 %1:numel(keys) %ceil(0.75*numel(keys)):numel(keys)
=======
for fnum = 1:numel(keys) %ceil(0.75*numel(keys)):numel(keys)
>>>>>>> Stashed changes
    matfile = [keys{fnum}(1:end-3) 'mat'];
    load(matfile);
    
    h=figure('Visible','off','Position',[77          91        2215         420]);
    subplot(10,1,2:10);
<<<<<<< Updated upstream
    imagesc(t,f(7:172),flipud(s(7:172,:)));
    %colormap(1-hot);
=======
    imagesc(t,f(7:172),s(7:172,:));
>>>>>>> Stashed changes
    hold on;
    lbls = nan*zeros(size(s));
    phrases = return_phrase_times(elements{fnum});
    for phrasenum = 1:numel(phrases.phraseType)
        tonset = phrases.phraseFileStartTimes(phrasenum);
        toffset = phrases.phraseFileEndTimes(phrasenum);
        line([tonset tonset],[freq_min freq_max],'Color',[1 1 1],'LineStyle','--');
        line([toffset toffset],[freq_min freq_max],'Color',[0.5 0.5 0.5],'LineStyle','--'); 
    end
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    subplot(10,1,1);
    
    for phrasenum = 1:numel(phrases.phraseType)
        tonset = phrases.phraseFileStartTimes(phrasenum);
        toffset = phrases.phraseFileEndTimes(phrasenum);
        plot(t((t>=tonset) & (t<=toffset)),ones(1,sum((t>=tonset) & (t<=toffset))),'Color', ...
            colors(find(syllables == phrases.phraseType(phrasenum)),:),'LineWidth',10);
        hold on;
    end
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    xlim([0 t(end)]);
    set(gca,'color','none');
    axis off;
    ylim([0.9 2]);
    tokens = regexp(matfile,'_','split');
    title(['bird: ' tokens{1} ', file: ' tokens{2}]);
    saveas(h,fullfile(targetdir,[keys{fnum}(1:end-3) 'png']));
    hgclose(h);
    display(fnum/numel(keys));
    display(matfile)
end    

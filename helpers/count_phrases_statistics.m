%% phrase number statistics
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
cd ('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
load('lrb85315auto_annotation3.mat');
numphrases = [];
for fnum = 1:numel(keys)
    phrases = return_phrase_times(elements{fnum});
    numphrases = [numphrases; numel(phrases.phraseType)];
end

%% plot
[n,x] = hist(numphrases,1:25);
figure;
plot(x,n/sum(n));
xlabel('# phrases'); ylabel('Prob.')
set(gca,'FontSize',24);
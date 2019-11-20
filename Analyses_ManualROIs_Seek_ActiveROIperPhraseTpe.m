%% Pre-requisits
% This should be run after running
% [resmat, state_count, state_labels] = create_first_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
% [trans2, state_labels] = create_second_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
% locate_roi_phrased_locked_activity_types (shoud be turned into a function
% ...

%% repositories
addpath(genpath('/Users/yardenc/Documents/GitHub/BirdSongBout'),'-end');
addpath(genpath('/Users/yardenc/Documents/GitHub/pmtk3'),'-end');
addpath(genpath('/Users/yardenc/Documents/GitHub/VideoAnalysisPipeline'),'-end');
GithubDir = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/';
CNMFEfolder = [GithubDir 'CNMF_E'];
addpath(genpath(CNMFEfolder),'-end');
%%
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009auto_annotation1_fix' 'baseROIdata_'};
birdnum = 3;
switch birdnum
    case 1
        bird_params = bird1_params;
    case 2
        bird_params = bird2_params;
    case 3
        bird_params = bird3_params;
end
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
file_prefix = bird_params{5}; 
%%
switch birdnum
    case 1
        % For lrb85315:
        ignore_dates = {'2017_04_19'};
        ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
        join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
        path_to_annotation_file = fullfile('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav','lrb85315auto_annotation5_fix.mat');
        to_normalize = 1;
        include_zero = 1;
        % Leads to 30 phrase types
    case 2
        % For lbr3022
        ignore_dates = {'2017_04_05' '2017_04_06' '2017_04_11' '2017_04_12' '2017_04_14' ...
        '2017_04_16' '2017_04_20' '2017_04_21' '2017_04_23' '2017_04_25' '2017_04_26' ...
        '2017_04_27' '2017_04_30' '2017_05_03'};
        ignore_entries = [-1 100 102 101 103];
        join_entries = {[2 206]}; %{[207 307 407] [404 405] [208 209] [200 309]};
        path_to_annotation_file = fullfile('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lbr3022/movs/wav','lbr3022auto_annotation5_alexa.mat');
        to_normalize = 1;
        include_zero = 1;
    case 3
        % For lbr3009
        ignore_dates = {'2017_04_27' '2017_04_28' '2017_06_05' '2017_06_06' '2017_06_07' '2017_06_08' '2017_06_09' '2017_06_12' '2017_06_13' ... 
            '2017_06_14' '2017_06_15' '2017_06_16' '2017_06_19' '2017_06_20' '2017_06_21' '2017_06_22' '2017_06_27' ...
            '2017_06_28' '2017_06_29' '2017_06_30' '2017_07_03' '2017_07_04' '2017_07_06' '2017_07_07' '2017_07_10' ...
            '2017_07_11' '2017_07_12' '2017_07_13' '2017_07_14' '2017_07_18' '2017_07_19' '2017_07_20' '2017_07_21'};
        ignore_entries = [-1 100 102 101 103];
        join_entries = {}; %{[207 307 407] [404 405] [208 209] [200 309]};
        path_to_annotation_file = fullfile('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lbr3009/movs/wav','lbr3009auto_annotation1_fix.mat');
        to_normalize = 1;
        include_zero = 0;
end





%%
clc;
 [resmat, state_count, state_labels] = create_first_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero);
 [trans2, ~] = create_second_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero);
 [results,syllables] = locate_roi_phrased_locked_activity_types(birdnum,ignore_dates,ignore_entries,join_entries,include_zero,'spikes',1);
 disp('done');
 
 
 %% check in which phrase type the ROI is active
 %% create directories
cd(['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/']);
required_directories = {'SignalInPhraseTypes'};
for dirnum = 1:numel(required_directories)
    if ~exist(required_directories{dirnum},'dir')
        mkdir(required_directories{dirnum});
    end
end

%% 1-1 use the estimated spike rates 's'
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/SignalInPhraseTypes'];
signal_thr = 0.05;
hmm_thr = 0.005;
%clc;
display('Signals in different phrases:');
activity_in_phrases = [];
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        tmp = zeros(1,numel(syllables)); 
        vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*...
            (results(daynum).Max(roinum,:) >= results(daynum).Max_before(roinum,:)).*...
            (results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        tmp(vec) = results(daynum).hmm(roinum,vec); %1
%         (results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*...
%                    (results(daynum).Max(roinum,:) >= results(daynum).Max_before(roinum,:)).*...
        activity_in_phrases = [activity_in_phrases; tmp];
            % previous phrase type
            
           
            
      
    end
end
%% calculate ** specificity ** distribution
x = nansum(activity_in_phrases');
y = activity_in_phrases'./(ones(size(activity_in_phrases',1),1)*x);
thr = 0.05;
thr2 = 0.9;
y = y(:,x>thr);
res = zeros(1,size(y,2));
for i = 1:numel(res)
    z = y(:,i);
    z = sort(z,'descend');
    res(i) = min(find(cumsum(z)>thr2));
end
%% Plot

%%
cd('/Users/yardenc/Documents/Projects/CohenGardner2017_CanaryHVCImaging/Code');
load Spike_based_syllables_per_roi;
figure;
plot(0:10,mean(res25)); hold on;
plot(0:10,mean(res50)); 
plot(0:10,mean(res75)); 
set(gca,'YScale','log')

mn25 = mean(res25); mn25(1) = 0; mn25 = mn25/sum(mn25);
mn50 = mean(res50); mn50(1) = 0; mn50 = mn50/sum(mn50);
mn75 = mean(res75); mn75(1) = 0; mn75 = mn75/sum(mn75);

figure;
%semilogy(0:10,[mn25/mn25(2);mn50/mn50(2);mn75/mn75(2)]','o'); hold on;
semilogy(0:10,[mn25;mn50;mn75]','o');
% bar(0:10,mn50); 
% bar(0:10,mn75); 
%%
cd('/Users/yardenc/Documents/Projects/CohenGardner2017_CanaryHVCImaging/Code');
load Spike_based_syllables_per_roi_counts;
t25 = sum(n25); p25 = [0 t25(2:end)]/sum(t25(2:end)); 
t50 = sum(n50); p50 = [0 t50(2:end)]/sum(t50(2:end)); 
t75 = sum(n75); p75 = [0 t75(2:end)]/sum(t75(2:end)); 

% shrink the y axis for graphical reasons
%p25(2) = p25(2)-0.4;
%p50(2) = p50(2)-0.4;
%p75(2) = p75(2)-0.4;

figure; 
plot(0:10,p25/p25(2),'Marker','.','MarkerSize',15,'LineStyle','none'); hold on;
plot(0:10,p50/p50(2),'Marker','.','MarkerSize',15,'LineStyle','none');
plot(0:10,p75/p75(2),'Marker','.','MarkerSize',15,'LineStyle','none');
ylim([0 0.3]); xlim([1.5 4.5]); xticks(2:4);
yticks([0.1 0.2]);
set(gca,'FontSize',32);
xlabel('# syllables');
ylabel('Ratio to single-locked')

%% an alternative: use the integral
 [results,syllables] = per_phrase_activity_mass(birdnum,ignore_dates,ignore_entries,join_entries,include_zero,'spikes',0);
%%
fraction_thr = 0.25;
%clc;
display('Signals reliability in different phrases by integral:');
activity_in_phrases = [];
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).res50,1)
        activity_in_phrases = [activity_in_phrases; sum(results(daynum).res50(roinum,:) > fraction_thr)];
    end
end
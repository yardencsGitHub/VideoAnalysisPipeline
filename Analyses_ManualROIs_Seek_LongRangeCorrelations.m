%% Pre-requisits
% This should be run after running
% [resmat, state_count, state_labels] = create_first_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
% [trans2, state_labels] = create_second_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
% locate_roi_phrased_locked_activity_types (shoud be turned into a function
% ...
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009auto_annotation1_fix' 'baseROIdata_'};
birdnum = 1;
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
 [resmat, state_count, state_labels] = create_first_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero);
 [trans2, ~] = create_second_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero);
 [results,syllables] = locate_roi_phrased_locked_activity_types(birdnum,ignore_dates,ignore_entries,join_entries,include_zero);
%% Analysis - seek activity-centric correlations
% 1. look for phrases in which there is activity (hmm + high max)
% edges should be a parameter,
% Another option is have it (edges) depend on the hmm
% 2. check corr to durations and identity of next and previous phrases -
% alert findings
% 3.  check high signal in prev\post phrase (high max) --> alert findings
% 4.  check incoming and outgoing phrase types --> alert transition types
% (converging, diverging, ballistic .. output the actual transitions -->
% list)
% 5.  for each transition type. If going over a certain occurance threshold
% (say 0.1) condition on transition and look for 2nd order correlations of
% durations (for each 2nd order type) and of phrase type.
% For each correlation .. create a figure of the curves and the statistics

%%
% [results, syllables] = locate_roi_phrased_locked_activity_types(1,[],[],[],[],[]);
% cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs');
% [resmat, state_count, state_labels] = create_first_order_transition_matrix('lrb85315auto_annotation5_fix.mat',{'2017_04_19'},[-1 100 102 101 103 202 406 408 409 402 403],{[207 307 407] [404 405] [208 209] [200 309]},1,1);
%%
% for lrb85315
% ignore_dates = {'2017_04_19'};
% 
% ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
% 
% join_entries = {[207 307 407] [404 405] [208 209] [200 309]};

%% create directories
cd(['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations']);
required_directories = {'PostIdVsSignal' 'PrevIdVsSignal' 'Post2ndIdVsSignal' 'Prev2ndIdVsSignal' ...
    'PostDurationVsSignal' 'PrevDurationVsSignal' 'Post2ndDurationVsSignal' 'Prev2ndDurationVsSignal' ...
    'Sustained'};
for dirnum = 1:numel(required_directories)
    if ~exist(required_directories{dirnum},'dir')
        mkdir(required_directories{dirnum});
    end
end

%% 1-1 look for phrases in which there is activity (hmm + high max) and correlated to post phrase identity
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/PostIdVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to next phrase:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            
            h=figure('Position',[360    61   560   637],'Visible','off');
            ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[nan sylnum nan],roinum,1,0,[-3], ...
                'edges',[0 0],'use_residuals',0,'display_opt',0,'bird_number',3);
            if p<0.05
%                 f = gcf;
%                 h1 = get(ax,'Parent');
%                 fx = get(f,'Children');
                figure(h);
                h2 = subplot(2,1,2);
                c = copyobj(get(hndls(2),'Children'),h2);
                title(['1 way ANOVA: F =' num2str(r) ', p = ' num2str(p)]);
                xticks(1:numel(gnames)); xticklabels(gnames);
                display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                    num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                outputfilename = fullfile(targetdir,['signal2postId_anova_Date_' results(daynum).Date ...
                    '-roi' num2str(roinum) '-syllable' ...
                    num2str(sylnum) '.png']);
                saveas(h,outputfilename);
            end
            close all;
            
        end
    end
end

%% 1-2 look for phrases in which there is activity (hmm + high max) and correlated to prev phrase identity
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/PrevIdVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to previous phrase:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_before(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            
            h=figure('Position',[360    61   560   637],'Visible','off');
            ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[nan sylnum nan],roinum,1,0,[-1], ...
                'edges',[0 0],'use_residuals',0,'display_opt',0,'bird_number',3);
            if p<0.05
%                 f = gcf;
%                 h1 = get(ax,'Parent');
%                 fx = get(f,'Children');
                figure(h);
                h2 = subplot(2,1,2);
                c = copyobj(get(hndls(2),'Children'),h2);
                title(['1 way ANOVA: F =' num2str(r) ', p = ' num2str(p)]);
                xticks(1:numel(gnames)); xticklabels(gnames);
                display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                    num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                outputfilename = fullfile(targetdir,['signal2prevId_anova_Date_' results(daynum).Date ...
                    '-roi' num2str(roinum) '-syllable' ...
                    num2str(sylnum) '.png']);
                saveas(h,outputfilename);
            end
            close all;
            
        end
    end
end
%% 2-1 look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) and correlated to
% post phrase durations
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/PostDurationVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to post phrase duration:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            post_types = find(resmat(state_labels == sylnum,:) >= 0.05);
            for post_cnt = 1:numel(post_types)
                post_sylnum = state_labels(post_types(post_cnt));
                h=figure('Position',[360    61   560   637],'Visible','off');
                ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[nan sylnum post_sylnum],roinum,1,0,[3], ...
                'edges',[2 0],'use_residuals',0,'display_opt',0,'bird_number',3);
                if p<0.05
%                     f = gnames;
%                     h1 = get(ax,'Parent');
%                     fx = get(f,'Children');
%                     figure(h1);
%                     h2 = subplot(2,1,2);
%                     c = copyobj(get(fx,'Children'),h2);
                    figure(h);
                    h2 = subplot(2,1,2);
                    c = copyobj(get(hndls(2),'Children'),h2);
                    set(h2,'FontSize',16);
                    title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);
                    
                    display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                        num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                    outputfilename = fullfile(targetdir,['signal2post_duration_pearson_Date_' results(daynum).Date ...
                        '-roi' num2str(roinum) '-syllables ' ...
                        num2str(sylnum) '-->' num2str(post_sylnum) '.png']);
                    saveas(h,outputfilename);
                end
                close all;

            end
        end
    end
end
       
%% 2-2 look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) and correlated to
% prev phrase durations
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/PrevDurationVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to prev phrase duration:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_before(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            prev_types = find(resmat(:,state_labels == sylnum) >= 0.1);
            for prev_cnt = 1:numel(prev_types)
                prev_sylnum = state_labels(prev_types(prev_cnt));
                h=figure('Position',[360    61   560   637],'Visible','off');
                ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[prev_sylnum sylnum nan],roinum,1,0,[1], ...
                'edges',[0 2],'use_residuals',0,'display_opt',0,'bird_number',3);
                if p<0.05
%                     f = gnames;
%                     h1 = get(ax,'Parent');
%                     fx = get(f,'Children');
%                     figure(h1);
%                     h2 = subplot(2,1,2);
%                     c = copyobj(get(fx,'Children'),h2);
                    figure(h);
                    h2 = subplot(2,1,2);
                    c = copyobj(get(hndls(2),'Children'),h2);
                    set(h2,'FontSize',16);
                    title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);
                    
                    display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                        num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                    outputfilename = fullfile(targetdir,['signal2rev_duration_pearson_Date_' results(daynum).Date ...
                        '-roi' num2str(roinum) '-syllables ' ...
                        num2str(prev_sylnum) '-->' num2str(sylnum) '.png']);
                    saveas(h,outputfilename);
                end
                close all;

            end
        end
    end
end
%% 3-1 look for phrases in which there is activity (hmm + high max) 
% % AND go over transition types (with occurance > 0.1) and correlated to
% % 2n post phrase identity
trimmed_labels = state_labels(2:end-1);
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/Post2ndIdVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to 2nd post phrase:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            post_types = find(resmat(state_labels == sylnum,:) >= 0.05);
            post_types = setdiff(post_types,[1 numel(state_labels)]);
            for post_cnt = 1:numel(post_types)
                post_sylnum = state_labels(post_types(post_cnt));
                h=figure('Position',[360    61   560   637],'Visible','off');
                ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[nan sylnum post_sylnum nan],roinum,1,0,[-4], ...
                    'edges',[0 0.2],'use_residuals',0,'display_opt',0,'bird_number',3);
                if p<0.05
                    xlim(ax,[-2 5]);
%                     f = gcf;
%                     h1 = get(ax,'Parent');
%                     fx = get(f,'Children');
%                     figure(h1);
%                     h2 = subplot(2,1,2);
%                     c = copyobj(get(fx,'Children'),h2);
                    figure(h);
                    h2 = subplot(2,1,2);
                    c = copyobj(get(hndls(2),'Children'),h2);
                    title(['1 way ANOVA: F =' num2str(r) ', p = ' num2str(p)]);
                    xticks(1:numel(gnames)); xticklabels(gnames); 
                    display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                        num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                    outputfilename = fullfile(targetdir,['signal2post_2nd_Id_anova_Date_' results(daynum).Date ...
                        '-roi' num2str(roinum) '-syllable' ...
                        num2str(sylnum) '-->' num2str(post_sylnum) '.png']);

                    saveas(h,outputfilename);
                end
                close all;

            end
        end
    end
end

%% 3-2 look for phrases in which there is activity (hmm + high max) 
% % AND go over transition types (with occurance > 0.1) and correlated to
% % 2n prev phrase identity
trimmed_labels = state_labels(2:end-1);
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/Prev2ndIdVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to 2nd prev phrase:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            prev_types = find(resmat(:,state_labels == sylnum) >= 0.05);
            prev_types = setdiff(prev_types,[1 numel(state_labels)]);
            for prev_cnt = 1:numel(prev_types)
                prev_sylnum = state_labels(prev_types(prev_cnt));
                h=figure('Position',[360    61   560   637],'Visible','off');
                ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,3,[nan prev_sylnum sylnum nan],roinum,1,0,[-1], ...
                    'edges',[0.2 0],'use_residuals',0,'display_opt',0,'bird_number',3);
                if p<0.05
                    xlim(ax,[-5 2]);
%                     f = gcf;
%                     h1 = get(ax,'Parent');
%                     fx = get(f,'Children');
%                     figure(h1);
%                     h2 = subplot(2,1,2);
%                     c = copyobj(get(fx,'Children'),h2);
                    figure(h);
                    h2 = subplot(2,1,2);
                    c = copyobj(get(hndls(2),'Children'),h2);
                    title(['1 way ANOVA: F =' num2str(r) ', p = ' num2str(p)]);
                    xticks(1:numel(gnames)); xticklabels(gnames); 
                    display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                        num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                    outputfilename = fullfile(targetdir,['signal2post_2nd_Id_anova_Date_' results(daynum).Date ...
                        '-roi' num2str(roinum) '-syllable' ...
                        num2str(prev_sylnum) '-->' num2str(sylnum) '.png']);

                    saveas(h,outputfilename);
                end
                close all;

            end
        end
    end
end
%% 4-1 look for phrases in which there is activity (hmm + high max) 
%AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS and correlated to
%prev phrase durations. This is done with the 2nd order transition probabilities 
%targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PostDurationVsSignal';
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/PrevDurationVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to previous phrase duration:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr)); %(results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            prev_types = find(resmat(:,state_labels == sylnum) >= 0.05);   
            %prev_types = setdiff(prev_types,1);
            for prev_cnt = 1:numel(prev_types)
                prev_sylnum = state_labels(prev_types(prev_cnt));
                post_types = find(trans2(state_labels == prev_sylnum,state_labels == sylnum,:) >= 0.05);
                
                for post_cnt = 1:numel(post_types)
                    post_sylnum = state_labels(post_types(post_cnt));
                    h=figure('Position',[360    61   560   637],'Visible','off');
                    ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[prev_sylnum sylnum post_sylnum],roinum,1,0,[1],...
                        'edges',[0 2],'use_residuals',0,'display_opt',0,'bird_number',3);
                    if p<0.05
                        figure(h);
                        h2 = subplot(2,1,2);
                        c = copyobj(get(hndls(2),'Children'),h2);
                        set(h2,'FontSize',16);
                        title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);

                        display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                            num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                        outputfilename = fullfile(targetdir,['signal2prev_duration_pearson_Date_' results(daynum).Date ...
                            '-roi' num2str(roinum) '-syllables ' ...
                            num2str(prev_sylnum) '-->' num2str(sylnum) '-->' num2str(post_sylnum) '.png']);
                        saveas(h,outputfilename);
                    end
                    close all;
                end
            end
        end
    end
end

%% 4-2 look for phrases in which there is activity (hmm + high max) 
%AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS and correlated to
%post phrase durations. This is done with the 2nd order transition probabilities 
%targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PostDurationVsSignal';
targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
    '/ManualROIs/LongRangeCorrelations/PostDurationVsSignal'];
signal_thr = 0.05;
hmm_thr = 0.05;
clc;
display('Signal relation to next phrase duration:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr)); %(results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            prev_types = find(resmat(:,state_labels == sylnum) >= 0.05);   
            %prev_types = setdiff(prev_types,1);
            for prev_cnt = 1:numel(prev_types)
                prev_sylnum = state_labels(prev_types(prev_cnt));
                post_types = find(trans2(state_labels == prev_sylnum,state_labels == sylnum,:) >= 0.05);
                
                for post_cnt = 1:numel(post_types)
                    post_sylnum = state_labels(post_types(post_cnt));
                    h=figure('Position',[360    61   560   637],'Visible','off');
                    ax = subplot(2,1,1); [hndls,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[prev_sylnum sylnum post_sylnum],roinum,1,0,[3],...
                        'edges',[2 0],'use_residuals',0,'display_opt',0,'bird_number',3);
                    if p<0.05
                        figure(h);
                        h2 = subplot(2,1,2);
                        c = copyobj(get(hndls(2),'Children'),h2);
                        set(h2,'FontSize',16);
                        title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);

                        display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                            num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                        outputfilename = fullfile(targetdir,['signal2post_duration_pearson_Date_' results(daynum).Date ...
                            '-roi' num2str(roinum) '-syllables ' ...
                            num2str(prev_sylnum) '-->' num2str(sylnum) '-->' num2str(post_sylnum) '.png']);
                        saveas(h,outputfilename);
                    end
                    close all;
                end
            end
        end
    end
end



%% 5 STAT correction look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS and correlated to
% prev/post phrase durations. This is done with the 2nd order transition probabilities.
% USE RESIDUALS (change is made in
% LongRangeLockedSingleDayManualROIs_function.m
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PrevDurationVsSignal';
% targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
%      '/ManualROIs/LongRangeCorrelations/PostDurationVsSignal'];
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Signal relation to previous phrase:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
%         %(results(daynum).Max(roinum,:) <= results(daynum).Max_before(roinum,:)).*
%         for syl = 1:numel(vec)
%             sylnum = syllables(vec(syl));
%             % previous phrase type
%             prev_types = find(resmat(:,state_labels == sylnum) >= 0.1);   
%             %prev_types = setdiff(prev_types,1);
%             for prev_cnt = 1:numel(prev_types)
%                 prev_sylnum = state_labels(prev_types(prev_cnt));
%                 post_types = find(trans2(state_labels == prev_sylnum,state_labels == sylnum,:) >= 0.1);
%                 
%                 for post_cnt = 1:numel(post_types)
%                     post_sylnum = state_labels(post_types(post_cnt));
%                     h=figure('Position',[360    61   560   637],'Visible','off');
%                     ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[prev_sylnum sylnum post_sylnum],roinum,1,0,[3], ...
%                         'edges',[1 0],'use_residuals',1,'display_opt',0);
%                     if p<0.05
%                         f = gnames;
%                         h1 = get(ax,'Parent');
%                         fx = get(f,'Children');
%                         figure(h1);
%                         h2 = subplot(2,1,2);
%                         c = copyobj(get(fx,'Children'),h2);
%                         set(h2,'FontSize',16);
%                         title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);
% 
%                         display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                             num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
%                         outputfilename = fullfile(targetdir,['STAT_signal2post_duration_pearson_Date_' results(daynum).Date ...
%                             '-roi' num2str(roinum) '-syllables ' ...
%                             num2str(prev_sylnum) '-->' num2str(sylnum) '-->' num2str(post_sylnum) '.png']);
%                         saveas(h,outputfilename);
%                     end
%                     close all;
%                 end
%             end
%         end
%     end
% end
   
%% 6 STAT correction look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS and correlated to
% 2nd post phrase durations. This is done with the 2nd order transition probabilities.
% USE RESIDUALS (change is made in
% LongRangeLockedSingleDayManualROIs_function.m
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/Post2ndDurationVsSignal';
% targetdir = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name ...
%      '/ManualROIs/LongRangeCorrelations/Post2ndDurationVsSignal'];
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Signal relation to 2nd post phrase duration:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
%         %(results(daynum).Max(roinum,:) <= results(daynum).Max_before(roinum,:)).*
%         for syl = 1:numel(vec)
%             sylnum = syllables(vec(syl));
%             % previous phrase type
%             post1_types = find(resmat(state_labels == sylnum,:) >= 0.1);   
%             %prev_types = setdiff(prev_types,1);
%             for post1_cnt = 1:numel(post1_types)
%                 post1_sylnum = state_labels(post1_types(post1_cnt));
%                 post2_types = find(trans2(state_labels == sylnum,state_labels == post1_sylnum,:) >= 0.1);
%                 
%                 for post_cnt = 1:numel(post2_types)
%                     post2_sylnum = state_labels(post2_types(post_cnt));
%                     h=figure('Position',[360    61   560   637],'Visible','off');
%                     ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,1,[sylnum post1_sylnum post2_sylnum],roinum,1,0,[3], ...
%                     'edges',[1 0],'use_residuals',1,'display_opt',0);
%                     if p<0.05
%                         f = gnames;
%                         h1 = get(ax,'Parent');
%                         fx = get(f,'Children');
%                         figure(h1);
%                         h2 = subplot(2,1,2);
%                         c = copyobj(get(fx,'Children'),h2);
%                         set(h2,'FontSize',16);
%                         title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);
% 
%                         display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                             num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
%                         outputfilename = fullfile(targetdir,['STAT_signal2post2nd_duration_pearson_Date_' results(daynum).Date ...
%                             '-roi' num2str(roinum) '-syllables ' ...
%                             num2str(sylnum) '-->' num2str(post1_sylnum) '-->' num2str(post2_sylnum) '.png']);
%                         saveas(h,outputfilename);
%                     end
%                     close all;
%                 end
%             end
%         end
%     end
% end

%% 7 look for sustained activity - consecutive phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS . This is done with the 2nd order transition probabilities.
% cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs');
% load results2018_01_15; % contains the transition statistics and state, syllable variables
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/Sustained';
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Sustained Signals:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr).*...
%             (results(daynum).Max_before(roinum,:) >= signal_thr).*(results(daynum).Max_after(roinum,:) >= signal_thr));
%             for syl = 1:numel(vec)
%                 sylnum = syllables(vec(syl));
%                 h=figure('Position',[360    61   560   637],'Visible','off');
%                 ax = axes; [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,1,[nan sylnum nan],roinum,1,0,[2]);
%                 display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                             num2str(sylnum) ', sustained' ]);
%                 outputfilename = fullfile(targetdir,['Sustained_duration_pearson_Date_' results(daynum).Date ...
%                             '-roi' num2str(roinum) '-syllable ' ...
%                             num2str(sylnum) '.png']);
%                 saveas(h,outputfilename);
%             
%                 close all;
%             end
%         
%     end
% end
        
%% 8 sustained in a different way
% cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/ManualROIs');
% load results2018_01_15; % contains the transition statistics and state, syllable variables
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/Sustained';
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Sustained Signals:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
%         %(results(daynum).Max(roinum,:) <= results(daynum).Max_before(roinum,:)).*
%         for syl = 1:numel(vec)
%             sylnum = syllables(vec(syl));
%             % previous phrase type
%             post1_types = find(resmat(state_labels == sylnum,:) >= 0.1); 
%             post1_sylnums = state_labels(post1_types);
%             post1_sylnums = setdiff(post1_sylnums,[1000 -1000]);
%             post1_sylnums(results(daynum).Max(roinum,ismember(state_labels(2:end-1),post1_sylnums)) < 0.1) = [];
%             %prev_types = setdiff(prev_types,1);
%             for post1_cnt = 1:numel(post1_sylnums)
%                 post1_sylnum = post1_sylnums(post1_cnt);
%                 post2_types = find(trans2(state_labels == sylnum,state_labels == post1_sylnum,:) >= 0.1);
%                 post2_sylnums = state_labels(post2_types);
%                 post2_sylnums = setdiff(post2_sylnums,[1000 -1000]);
%                 post2_sylnums(results(daynum).Max(roinum,ismember(state_labels(2:end-1),post2_sylnums)) < 0.1) = [];
%                 %post2_types(results(daynum).Max(roinum,post2_types) < 0.1) = [];
%                 for post_cnt = 1:numel(post2_sylnums)
%                     post2_sylnum = post2_sylnums(post_cnt);
%                     h=figure('Position',[360    61   560   637],'Visible','off');
%                     ax = axes; [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,1,[sylnum post1_sylnum post2_sylnum],roinum,1,0,[1 2 3]);
%                     if (gnames >= 5)
%                         display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                                 num2str(sylnum) ', sustained' ]);
%                         outputfilename = fullfile(targetdir,['Sustained_duration_n2_Date_' results(daynum).Date ...
%                                 '-roi' num2str(roinum) '-syllables ' ...
%                                 num2str(sylnum) '-->' num2str(post1_sylnum) '-->' num2str(post2_sylnum) '.png']);
%                         saveas(h,outputfilename);
%                         close all;
%                     end
%                 end
%             end
%         end
%     end
% end

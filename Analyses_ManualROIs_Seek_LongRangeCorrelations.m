%% Pre-requisits
% This should be run after running
% [resmat, state_count, state_labels] = create_first_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
% [trans2, state_labels] = create_second_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
% locate_roi_phrased_locked_activity_types (shoud be turned into a function
% ...
% For lrb85315:
%   ignore_dates = {'2017_04_19'};
%   ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
%   join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
% Leads to 30 phrase types

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
ignore_dates = {'2017_04_19'};

ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];

join_entries = {[207 307 407] [404 405] [208 209] [200 309]};

%% 1 look for phrases in which there is activity (hmm + high max) and correlated to prev/post phrase identity
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PrevIdVsSignal';
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Signal relation to previous phrase:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_before(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
%         for syl = 1:numel(vec)
%             sylnum = syllables(vec(syl));
%             % previous phrase type
%             
%             h=figure('Position',[360    61   560   637],'Visible','off');
%             ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[nan sylnum nan],roinum,1,0,[-1]);
%             if p<0.05
%                 f = gcf;
%                 h1 = get(ax,'Parent');
%                 fx = get(f,'Children');
%                 figure(h1);
%                 h2 = subplot(2,1,2);
%                 c = copyobj(get(fx,'Children'),h2);
%                 title(['1 way ANOVA: F =' num2str(r) ', p = ' num2str(p)]);
%                 xticks(1:numel(gnames)); xticklabels(gnames);
%                 display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                     num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
%                 outputfilename = fullfile(targetdir,['signal2prevId_anova_Date_' results(daynum).Date ...
%                     '-roi' num2str(roinum) '-syllable' ...
%                     num2str(sylnum) '.png']);
%                 saveas(h,outputfilename);
%             end
%             close all;
%             
%         end
%     end
% end

%% 2 look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) and correlated to
% prev/post phrase durations
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PostDurationVsSignal';
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Signal relation to previous phrase:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
%         for syl = 1:numel(vec)
%             sylnum = syllables(vec(syl));
%             % previous phrase type
%             post_types = find(resmat(state_labels == sylnum,:) >= 0.1);
%             for post_cnt = 1:numel(post_types)
%                 post_sylnum = state_labels(post_types(post_cnt));
%                 h=figure('Position',[360    61   560   637],'Visible','off');
%                 ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[nan sylnum post_sylnum],roinum,1,0,[3]);
%                 if p<0.05
%                     f = gnames;
%                     h1 = get(ax,'Parent');
%                     fx = get(f,'Children');
%                     figure(h1);
%                     h2 = subplot(2,1,2);
%                     c = copyobj(get(fx,'Children'),h2);
%                     set(h2,'FontSize',16);
%                     title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);
%                     
%                     display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                         num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
%                     outputfilename = fullfile(targetdir,['signal2post_duration_pearson_Date_' results(daynum).Date ...
%                         '-roi' num2str(roinum) '-syllables ' ...
%                         num2str(sylnum) '-->' num2str(post_sylnum) '.png']);
%                     saveas(h,outputfilename);
%                 end
%                 close all;
% 
%             end
%         end
%     end
% end
%        

% %% 3 look for phrases in which there is activity (hmm + high max) 
% % AND go over transition types (with occurance > 0.1) and correlated to
% % 2n prev/post phrase identity
% trimmed_labels = state_labels(2:end-1);
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/Prev2ndIdVsSignal';
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Signal relation to previous phrase:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
%         for syl = 1:numel(vec)
%             sylnum = syllables(vec(syl));
%             % previous phrase type
%             prev_types = find(resmat(:,state_labels == sylnum) >= 0.1);
%             prev_types = setdiff(prev_types,[1 32]);
%             for prev_cnt = 1:numel(prev_types)
%                 prev_sylnum = state_labels(prev_types(prev_cnt));
%                 h=figure('Position',[360    61   560   637],'Visible','off');
%                 ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,3,[nan prev_sylnum sylnum nan],roinum,1,0,[-1]);
%                 if p<0.05
%                     xlim(ax,[-5 2]);
%                     f = gcf;
%                     h1 = get(ax,'Parent');
%                     fx = get(f,'Children');
%                     figure(h1);
%                     h2 = subplot(2,1,2);
%                     c = copyobj(get(fx,'Children'),h2);
%                     title(['1 way ANOVA: F =' num2str(r) ', p = ' num2str(p)]);
%                     xticks(1:numel(gnames)); xticklabels(gnames); 
%                     display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
%                         num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
%                     outputfilename = fullfile(targetdir,['signal2prev_2nd_Id_anova_Date_' results(daynum).Date ...
%                         '-roi' num2str(roinum) '-syllable' ...
%                         num2str(prev_sylnum) '-->' num2str(sylnum) '.png']);
% 
%                     saveas(h,outputfilename);
%                 end
%                 close all;
% 
%             end
%         end
%     end

%% 4 look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS and correlated to
% prev/post phrase durations. This is done with the 2nd order transition probabilities 
% targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PostDurationVsSignal';
% signal_thr = 0.1;
% hmm_thr = 0.1;
% clc;
% display('Signal relation to previous phrase:');
% for daynum = 1:numel(results)
%     for roinum = 1:size(results(daynum).Max,1)
%         vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_after(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
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
%                     ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[prev_sylnum sylnum post_sylnum],roinum,1,0,[3]);
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
%                         outputfilename = fullfile(targetdir,['signal2post_duration_pearson_Date_' results(daynum).Date ...
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


%% 5 STAT correction look for phrases in which there is activity (hmm + high max) 
% AND go over transition types (with occurance > 0.1) WITH IDENTICAL FLANKERS and correlated to
% prev/post phrase durations. This is done with the 2nd order transition probabilities.
% USE RESIDUALS (change is made in
% LongRangeLockedSingleDayManualROIs_function.m
targetdir = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/ManualROIs/LongRangeCorrelations/PrevDurationVsSignal';
signal_thr = 0.1;
hmm_thr = 0.1;
clc;
display('Signal relation to previous phrase:');
for daynum = 1:numel(results)
    for roinum = 1:size(results(daynum).Max,1)
        vec = find((results(daynum).Max(roinum,:) >= results(daynum).Max_before(roinum,:)).*(results(daynum).Max(roinum,:) >= signal_thr).*(results(daynum).hmm(roinum,:) >= hmm_thr));
        for syl = 1:numel(vec)
            sylnum = syllables(vec(syl));
            % previous phrase type
            prev_types = find(resmat(:,state_labels == sylnum) >= 0.1);   
            %prev_types = setdiff(prev_types,1);
            for prev_cnt = 1:numel(prev_types)
                prev_sylnum = state_labels(prev_types(prev_cnt));
                post_types = find(trans2(state_labels == prev_sylnum,state_labels == sylnum,:) >= 0.1);
                
                for post_cnt = 1:numel(post_types)
                    post_sylnum = state_labels(post_types(post_cnt));
                    h=figure('Position',[360    61   560   637],'Visible','off');
                    ax = subplot(2,1,1); [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,results(daynum).Date,ignore_entries,join_entries,2,[prev_sylnum sylnum post_sylnum],roinum,1,0,[1]);
                    if p<0.05
                        f = gnames;
                        h1 = get(ax,'Parent');
                        fx = get(f,'Children');
                        figure(h1);
                        h2 = subplot(2,1,2);
                        c = copyobj(get(fx,'Children'),h2);
                        set(h2,'FontSize',16);
                        title(['Pearson: r =' num2str(r) ', p = ' num2str(p)]);

                        display([results(daynum).Date ' - roi #' num2str(roinum) ', syllable ' ...
                            num2str(sylnum) ', ANOVA (F,p) = ' num2str(r) ',' num2str(p)]);
                        outputfilename = fullfile(targetdir,['STAT_signal2prev_duration_pearson_Date_' results(daynum).Date ...
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
       
        
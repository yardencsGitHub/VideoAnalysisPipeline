%% Create PSTs for canaries lrb85314, lbr3022, lbr3009
% This script creates the PSTs and figures
% To reduce noise some rare syllables are removed (the variable 'ignore_entries') and the criterion for
% removal is described for each bird (since the number of files changes
% between them)
% some dates are ignored (the 'ignore_dates' list) because of data
% corruption.
% Some entries are joined (the list 'join_entries') because they belong to
% the same syllable class and were separated due to user OCD.. :)
% 'include_zero' is 0/1 for including the label '0' (in birds 853 and
% 3022), 'min_phrases' sets the minimal # of phrases per song bout to
% include
%% repos. requirements: pst, BirdSongBout
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/pst'),'-end');
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/BirdSongBout/BirdSongBout'),'-end');
%% lrb85315
path_to_annotation = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/lrb85315auto_annotation5_fix.mat';
ignore_dates = {'2017_04_18' '2017_04_19' '2017_04_20'};
ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
join_entries = {[207 307 407] [404 405] [200 309] [208 209]};
include_zero = 1;
min_phrases = 2;
load('/Users/yardenc/Documents/Projects/CohenGardner2017_CanaryHVCImaging/Figures/lrb853_colors.mat');
%%


% [DATA, syls] = convert_annotation_to_pst(path_to_annotation,{'2017_04_18' '2017_04_19' '2017_04_20'},[-1 100 102 101 103 202 406 408 409 402 403],{[207 307 407] [404 405] [200 309] [208 209]},1,2);
% 
% [F_MAT ALPHABET N PI]=pst_build_trans_mat(DATA,5);
% TREE = pst_learn(F_MAT,ALPHABET,N,'L',5);

%%
syl_labels = mat2cell(syls',ones(numel(syls),1),1);
%%
for rootnode = 1:numel(TREE(2).string)
    display(['Main branch: Syllble #' num2str(syls(TREE(2).string{rootnode}))]);
    figure; h_pie=pie(double(TREE(2).g_sigma_s(:,rootnode)),zeros(1,numel(syls)),cellfun(@num2str,syl_labels,'UniformOutput',0));
    for cnt = 2:2:numel(h_pie)
        h_pie(cnt).FontSize = 16;
        h_pie(cnt-1).FaceColor = colors(color_tags == str2num(h_pie(cnt).String),:);
    end
    title(['Root: ' num2str(syls(TREE(2).string{rootnode}))]);
    explore_branch(TREE,syls,2,rootnode);
end
    
function explore_branch(TREE,syls,level_num,branch_num)
    load('/Users/yardenc/Documents/Projects/CohenGardner2017_CanaryHVCImaging/Figures/lrb853_colors.mat');
    syl_labels = mat2cell(syls',ones(numel(syls),1),1);
   if (level_num < numel(TREE) & ~isempty(TREE(level_num + 1).parent))
       branches = find(TREE(level_num + 1).parent(1,:) == branch_num);
       for next_branch_num = 1:numel(branches)
           %display([level_num next_branch_num])
           display(['level ' num2str(level_num) ' Branch: ' num2str(syls(TREE(level_num + 1).string{branches(next_branch_num)}))]);
           figure; h_pie=pie(double(TREE(level_num+1).g_sigma_s(:,branches(next_branch_num))),zeros(1,numel(syls)),cellfun(@num2str,syl_labels,'UniformOutput',0));
            for cnt = 2:2:numel(h_pie)
                h_pie(cnt).FontSize = 16;
                h_pie(cnt-1).FaceColor = colors(color_tags == str2num(h_pie(cnt).String),:);
            end
            title(['Branch: ' num2str(syls(TREE(level_num + 1).string{branches(next_branch_num)}))]);
           explore_branch(TREE,syls,level_num + 1,branches(next_branch_num));
       end
       
   end
end
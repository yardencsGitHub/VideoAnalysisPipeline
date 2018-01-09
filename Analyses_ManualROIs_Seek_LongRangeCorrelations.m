%% Pre-requisits
% This should be run after running
% [resmat, state_count, state_labels] = create_first_order_transition_matrix(path_to_annotation_file,ignore_dates,ignore_entries,join_entries,to_normalize,include_zero)
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
% - check corr to durations and identity of next and previous phrases -
% alert findings
% - check high signal in prev\post phrase (high max) --> alert findings
% - check incoming and outgoing phrase types --> alert transition types
% (converging, diverging, ballistic .. output the actual transitions -->
% list)
% - for each transition type. If going over a certain occurance threshold
% (say 0.1) condition on transition and look for 2nd order correlations of
% durations (for each 2nd order type) and of phrase type.
% For each correlation .. create a figure of the curves and the statistics


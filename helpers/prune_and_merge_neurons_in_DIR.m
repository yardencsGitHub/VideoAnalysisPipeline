DIR = pwd;
file_list = dir('neuron_*.mat');
for fnum = 851:numel(file_list)
    fname = file_list(fnum).name;
    display(fname);
    load(fname);
    [A,C,C_raw,S] = prune_and_merge_neurons_postCNMFE(neuron);
    save(['pruned_' fname],'A','C','C_raw','S','-v7.3');
    close all;
end
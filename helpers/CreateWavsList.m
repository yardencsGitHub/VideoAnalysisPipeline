%% create the list of mov files to convert to matlab matrices
%cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
keysinfile = 'wavs_list';
%startfrom = 6427;
FILES = dir('*.wav');
ord = [];
for i = 1:numel(FILES)
    tokens = regexp(FILES(i).name,'_','split');
    ord = [ord; str2num(tokens{2})];
end
[locs,indx] = sort(ord);
startloc = find(locs == startfrom);
FILES = FILES(indx);
%%
wavs = {FILES(startloc:end).name};
save(keysinfile,'wavs');
%% create the list of mov files to convert to matlab matrices
cd('/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs');
keysinfile = 'FS_movies_list';
startfrom = 6427;
FILES = dir('*.mov');
ord = [];
for i = 1:numel(FILES)
    tokens = regexp(FILES(i).name,'_','split');
    ord = [ord; str2num(tokens{2})];
end
[locs,indx] = sort(ord);
startloc = find(locs == startfrom);
FILES = FILES(indx);
%%
keys = {FILES(startloc:end).name};
save(keysinfile,'keys');
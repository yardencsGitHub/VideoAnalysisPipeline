%%
cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
load lrb85315auto_annotation3;
%%
dates=[];
for fnum = 1:numel(keys);
    tokens = regexp(keys{fnum},'_','split');
    dates = [dates; datenum([tokens{3} '-' tokens{4} '-' tokens{5}])];
    
end

%%
datelist = unique(dates);
nums=[];
for i=1:numel(datelist);
    nums = [nums; sum(dates==datelist(i))];
end
%%
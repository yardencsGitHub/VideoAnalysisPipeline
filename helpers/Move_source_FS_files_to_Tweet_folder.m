%% move movies from storage folder, give them tweet names and check for duplicities
% assumes FS names in source and Tweet names in target

% last_idx = 7982;
% last_date = '2017_06_29';
% last_time = '08_24_46';
% bird_name = 'lrb85315';
%DamagedFolder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/too_large_or_damaged/';

DamagedMovs = dir([DamagedFolder ,'*.mov']);
TweetDateTime = {};
for fnum = 1: numel(DamagedMovs)
    tokens = regexp(DamagedMovs(fnum).name,'_','split');
    if strcmp(tokens{1},bird_name)   
        TweetDateTime = {TweetDateTime{:} DamagedMovs(fnum).name(end-22:end-4)};
    end
end
%%
%SourceFolder = '/Volumes/CanaryData/DATA/lrb853_15/movs/';
SourceMovs = dir([SourceFolder '*.mov']);
%TargetFolder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/';
TargetMovs = dir([TargetFolder ,'*.mov']);
TweetNums = [];
%TweetDateTime = {};
for fnum = 1: numel(TargetMovs)
    tokens = regexp(TargetMovs(fnum).name,'_','split');
    if strcmp(tokens{1},bird_name)
        TweetNums = [TweetNums; str2num(tokens{2})];
        TweetDateTime = {TweetDateTime{:} TargetMovs(fnum).name(end-22:end-4)};
    end
end


%%
cnt = 0;
for fnum = 1:numel(SourceMovs)
    sourcename = SourceMovs(fnum).name;
    try
        datenum(sourcename(1:10));
        tokens1 = regexp(sourcename(1:end-4),' ','split');
        tokens = regexp(tokens1{1},'-','split');

        datetimestr = char(join({tokens{:} tokens1{2:4}},'_'));
        flag = sum(cellfun(@(s)strcmp(s,datetimestr),TweetDateTime));
        if flag >= 1
            display(sourcename);
            cnt = cnt+1;
            display(cnt);
        else
           targetname = [ConvertFilenameFS2Tweet(sourcename,bird_name,last_idx+1) '.mov'];
           copyfile(fullfile(SourceFolder,sourcename),fullfile(TargetFolder,targetname));
           last_idx = last_idx+1;
        end
    catch em
        display('nodate');
        display(sourcename);
    end
end
    
%%
% DIR = pwd;
% load('lrb85315auto_annotation4.mat');
% 
% for fnum = 1:numel(keys)
%     t = regexp(keys{fnum},'_','split');
%     t{2} = sprintf('%04d',str2num(t{2}));
%     fname = char(join(t,'_'));
%     [SUCCESS,MESSAGE,MESSAGEID] = movefile(keys{fnum},fname);
%     keys{fnum} = fname;
%     elements{fnum}.filenum = t{2};
%     display(fname);
% end

%%
DIR = pwd;
% cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
% load('lrb85315auto_annotation4.mat');
cd(DIR);
d = dir('*.wav');
keys = {d.name};
for fnum = 1:numel(keys)
    try
        t = regexp(keys{fnum},'_','split');
        if numel(t{2})<4
            t{2} = sprintf('%04d',str2num(t{2}));
            fname = char(join(t,'_'));
            [SUCCESS,MESSAGE,MESSAGEID] = movefile(keys{fnum},fname);
            %keys{fnum} = fname;
            %elements{fnum}.filenum = t{2};
            display(fname);
        end
    catch em
    end
end
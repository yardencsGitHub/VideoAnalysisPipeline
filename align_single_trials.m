files = dir('SingleTrial*.mat');
for fnum = 1:1 %numel(files)
    load(files(fnum).name);
    [Y,d]  = AlignSingleTrial(Y);
    outname = ['aligned_' files(fnum).name];
    save(outname,'Y','Ysiz','-v7.3');
end
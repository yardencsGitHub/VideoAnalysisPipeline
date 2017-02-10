function AlignBatch(DIR,input_prefix,output_prefix,ref_mov_idx)
if (isempty(ref_mov_idx) || (ref_mov_idx == 0))
    ref_mov_idx = 1;
end

[FileName,PathName,FilterIndex] = uigetfile(fullfile(DIR,[input_prefix '*.mat']),'MultiSelect','on');
load(fullfile(PathName,FileName{ref_mov_idx}));
max_proj_ref = max(double(Y),[],3);
max_proj_ref = max_proj_ref - imgaussfilt(max_proj_ref,15);
%max_proj_ref = max_proj_ref-min(max_proj_ref(:));
%max_proj_ref = max_proj_ref/max(max_proj_ref(:));

[optimizer,metric] = imregconfig('monomodal');

for fnum = 1:numel(FileName)
    disp(fnum)
    load(fullfile(PathName,FileName{fnum}));
    max_proj_mov = max(double(Y),[],3);
    max_proj_mov = max_proj_mov - imgaussfilt(max_proj_mov,15);
    %max_proj_mov = max_proj_mov-min(max_proj_mov(:));
    %max_proj_mov = max_proj_mov/max(max_proj_mov(:));
    tform = imregtform(max_proj_mov,max_proj_ref,'rigid',optimizer,metric);
    d = tform.T(3,1:2);
    
    %if (max(abs(d))) < 20
        Y = imwarp(double(Y),tform,'OutputView',imref2d(size(max_proj_ref)));
    %else
    %    '*';
    %end
    outputfile = [output_prefix FileName{fnum}(length(input_prefix)+1:end)];
    save(outputfile,'Y','Ysiz','-v7.3');
    
end
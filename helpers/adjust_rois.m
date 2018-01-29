% This is a list of scripts that allow manipulating existing ROI files and
% ROI data.
% 1. Adjusting all roi shapes
% 2. Remove overlapping rois
% 3. create an image of a specific roi map (for examples)
% 4. Delete df/f files with a defined prefix (to get rid of old attempts)

cd('/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs/ROIdata');
basedir = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs/ROIdata';
dirs = dir;
dirs = dirs(4:end);
for dirnum = 26:numel(dirs)
    cd(fullfile(basedir,dirs(dirnum).name));
    filename = ['ROI_' dirs(dirnum).name];
    load(filename);
    f=figure; imagesc(ROI.reference_image); hx = gca;
    [nrows, ncolumns] = size(ROI.reference_image);
    [xi,yi]=meshgrid(1:ncolumns,1:nrows);
    colormap(gray);
    hold on;
    set(hx,'XTickLabel',[]);
    set(hx,'YTickLabel',[]);

    for el_num = 1:numel(ROI.coordinates)

       k=convhull(ROI.coordinates{el_num}(:,1),ROI.coordinates{el_num}(:,2));
       [aa]=ROI.coordinates{el_num}(k,:);
       plot(hx,aa(:,1),aa(:,2)); hold on;

    end
    for el_num = 1:numel(ROI.coordinates)
       k=convhull(ROI.coordinates{el_num}(:,1),ROI.coordinates{el_num}(:,2));
       [aa]=ROI.coordinates{el_num}(k,:);
       if size(aa,1)>10
           xorg = linspace(1,10,size(aa,1)-1);
           aa = interp1(xorg,aa(1:end-1,:),1:10);
       end
       h = impoly(hx,aa,'Closed',1);
       h1 = wait(h);
       delete(h);
       roi=inpolygon(xi,yi,h1(:,1),h1(:,2));
       [idx]=find(roi);
       ROI.coordinates{el_num} = [xi(idx) yi(idx)];
    end
    outname = ['new' filename];
    save(outname,'ROI');
    hgclose(f);
end
%% 2. now remove overlapping rois
cd('/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs/ROIdata');
basedir = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs/ROIdata';
dirs = dir;
dirs = dirs(4:end);
for dirnum = 1:numel(dirs)
    cd(fullfile(basedir,dirs(dirnum).name));
    filename = ['newROI_' dirs(dirnum).name '.mat'];
    if exist(filename)
        load(filename);
        tempROI.type = ROI.type;
        tempROI.reference_image = ROI.reference_image;
        tempROI.coordinates = {};
        tempROI.stats = [];
        mask = zeros(nrows,ncolumns);
        for roi_n = 1:numel(ROI.coordinates)
            roimask = zeros(480,640);
            roimask(ROI.coordinates{roi_n}(:,1)*nrows+ROI.coordinates{roi_n}(:,2))=1;
            if sum(sum(mask.*roimask)) == 0
                mask(ROI.coordinates{roi_n}(:,1)*nrows+ROI.coordinates{roi_n}(:,2))=1;
                tempROI.coordinates = {tempROI.coordinates{:} ROI.coordinates{roi_n}};
                tempROI.stats = [tempROI.stats, ROI.stats(roi_n)];
            end
        end
        display(filename); display([numel(ROI.coordinates) numel(tempROI.coordinates)]);
        ROI = tempROI;
        outname = ['nonoverlap_' filename];
        save(outname,'ROI');
    end
end
%% 3. create an image of a specific roi map (for examples)
mask = zeros(480,640);
nrows = 480;
textlocs = [];
for roi_n = 1:numel(ROI.coordinates)
   
    mask(ROI.coordinates{roi_n}(:,1)*nrows+ROI.coordinates{roi_n}(:,2))=roi_n;
    textlocs = [textlocs; mean(ROI.coordinates{roi_n}(:,1)) mean(ROI.coordinates{roi_n}(:,2))];
end
figure; imagesc(mask); hold on;
for roi_n = 1:numel(ROI.coordinates)
    text(textlocs(roi_n,1),textlocs(roi_n,2),num2str(roi_n),'FontSize',16,'HorizontalAlignment','center');
end
%% 4. Delete df/f files with a defined prefix (to get rid of old attempts)
prefix = 'testROIdata';
cd('/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs/ROIdata');
basedir = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/ManualROIs/ROIdata';
dirs = dir;
dirs = dirs(4:end);
for dirnum = 1:numel(dirs)
    cd(fullfile(basedir,dirs(dirnum).name));
    delete([prefix '*.mat']);
end

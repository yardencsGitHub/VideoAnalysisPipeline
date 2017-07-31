%%
syl = 11;
numrois = 20;
for roi = 1:numrois
    figure;
    for tr = 1:numel(rasters(syl,roi).data)
        plot(flipud(rasters(syl,roi).data{tr}),'k','LineWidth',2);
        hold on;
    end
    title(['ROI #' num2str(roi) ' ,syllable ' num2str(syl)]);
end
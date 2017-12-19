%%
threshold = 0.065;
syllable = 5;
flanking = [4 nan];
Day = 44;
cellnum = 16;
spikes = 2;
lockonset = 1;
qntthr = 0.999;

h = figure('Position',[588         234        1898        1104],'Visible','on');
subplot(3,4,1); plot(results(Day).Max(cellnum,:),'ro','MarkerSize',14,'MarkerFaceColor','r'); set(gca,'XTick',1:41); set(gca,'XTickLabel',syllables);
xlabel('syllables'); ylabel('0.9 quantile df/f');
xtickangle(90);
title(num2str(syltypes{Day}{cellnum}));
ax = subplot(3,4,[2 3 4 6 7 8]); [ax,r,p] = SequenceLockedSingleDayManualROIs_function(ax,results(Day).Date,syllable,flanking,cellnum,lockonset,spikes,0,2);
set(gca,'CameraPosition',[-14.5203 -192.6996    0.9790]);
[xidx,mn,se,mn2,se2,I1,durations] = MeanZscoredDFF_function(results(Day).Date,syllable,flanking,cellnum,lockonset,0,0,spikes);
subplot(3,4,11); plot(xidx,nanmean(I1));
xlabel('Real Time (Sec)'); ylabel('<df/f>')
activity_max = results(Day).Max(cellnum,syllables==syllable);
activity_min = results(Day).Min(cellnum,syllables==syllable);
hold on; line([0 0],[activity_min activity_max],'Color','r','LineStyle','--');
axis tight;

figure;
plot(xidx,nanmean(I1),'LineWidth',2);
xlabel('Real Time (Sec)','Color','k'); ylabel('<df/f>','Color','k');
hold on; line([0 0],[activity_min activity_max],'Color','r','LineStyle','--','LineWidth',2);
axis tight;
set(gca,'FontSize',32);

figure(h);
edges = [1 2];
fwhm = [];
cnts = [];
for cnt = 1:numel(durations)
    indx = find(xidx >=0-edges(1) & xidx <= durations(cnt)+edges(2));
    tmp = intersect(find(I1(cnt,:) > quantile(I1(cnt,xidx >=0 & xidx <= durations(cnt)),qntthr)/2),indx); %max(I1(cnt,xidx >=0 & xidx <= durations(cnt)))/2
    display([cnt xidx([min(tmp) max(tmp)])]);
    if isempty(tmp)
        fwhm = [fwhm; nan nan];
    else
        fwhm = [fwhm; xidx([min(tmp) max(tmp)])];
        if (fwhm(end,1) < 0)
            fwhm(end,1) = xidx(1+max(indx(I1(cnt,indx) < quantile(I1(cnt,xidx >=0 & xidx <= durations(cnt)),qntthr)/2 & ...
                xidx(indx) < 0)));
        end
        if (fwhm(end,2) > durations(cnt))
            fwhm(end,2) = xidx(min(indx(I1(cnt,indx) < quantile(I1(cnt,xidx >=0 & xidx <= durations(cnt)),qntthr)/2 & ...
                xidx(indx) > durations(cnt)))-1);
        end
    end
end

maxs = [];
for cnt = 1:numel(durations)
    indx = find(xidx >=0 & xidx <= durations(cnt));

    maxs = [maxs; max(I1(cnt,indx))];

end
subplot(3,4,5); plot(fwhm(:,:)); hold on; plot(durations(:));
xlabel('repetition #'); ylabel('Time from onset (sec)'); legend({'FWHM onset' 'FWHM offset' 'durations'});

figure;
plot(fwhm(:,:),'LineWidth',2); hold on; plot(durations(:),'LineWidth',2);
xlabel('repetition #','Color','k'); ylabel('Time from onset (sec)','Color','k'); legend({'FWHM onset' 'FWHM offset' 'durations'});
set(gca,'FontSize',32);
xlim([0 numel(durations)+1]);
figure(h);
repidx = find(fwhm(:,1)>-1.5 & maxs >= min(activity_max/2,0.1) & durations >= 0.1);
p = nan; r = p;
try
    [r,p] = corr([fwhm(repidx,:) diff(fwhm(repidx,:)')'],durations(repidx)); 
catch em
end
subplot(3,4,9); plot(maxs,'o','MarkerSize',14,'MarkerFaceColor','b'); 
hold on; plot(setdiff(1:numel(maxs),repidx),maxs(setdiff(1:numel(maxs),repidx)),'o','MarkerSize',14,'MarkerFaceColor',[0.5 0.5 0.5]); 
line([1 numel(maxs)],[activity_max/2.5 activity_max/2.5],'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','--')
xlabel('repetition #'); ylabel('df/f'); title(['Maxs ' num2str(r') ', ' num2str(p')]);

subplot(3,4,10); hist(maxs,20);
xlabel('Max value'); ylabel('# of repetitions');

fwhm = [];
cnts = [];
for cnt = 1:numel(durations)
    indx = find(xidx >=0-edges(1) & xidx <= durations(cnt)+edges(2));
    tmp = intersect(find(I1(cnt,:) > threshold),indx); %max(I1(cnt,xidx >=0 & xidx <= durations(cnt)))/2
    tmp1 = intersect(find(I1(cnt,:) <= threshold),indx);
    display([cnt xidx([min(tmp) max(tmp)])]);
    if isempty(tmp)
        fwhm = [fwhm; nan nan];
    else
        fwhm = [fwhm; xidx([min(tmp) max(tmp)])];
        try
        %fwhm = [fwhm; xidx([max(tmp1(xidx(tmp1)<durations(cnt))) max(tmp)])];
        catch em
           '-'; 
        end
        if (fwhm(end,1) < 0)
            try
            fwhm(end,1) = xidx(max(indx(I1(cnt,indx) < threshold & ...
                xidx(indx) < 0))+1);
            catch em
                '-';
            end
        end
        if (fwhm(end,2) > durations(cnt))
            fwhm(end,2) = xidx(min(indx(I1(cnt,indx) < threshold & ...
                xidx(indx) > durations(cnt)))-1);
        end
    end
end

maxs = [];
for cnt = 1:numel(durations)
    indx = find(xidx >=0 & xidx <= durations(cnt));

    maxs = [maxs; max(I1(cnt,indx))];

end
subplot(3,4,12); plot(fwhm(:,:)); hold on; plot(durations(:));
xlabel('repetition #'); ylabel('Time from onset (sec)'); legend({'Thr onset' 'Thr offset' 'durations'});
repidx = find(fwhm(:,1)>-1.5 & maxs >= min(activity_max/2,threshold) & durations >= 0.1);
p = nan; r = p;
try
    [r,p] = corr([fwhm(repidx,:) diff(fwhm(repidx,:)')'],durations(repidx)); 
catch em
end
title(['Maxs ' num2str(r') ', ' num2str(p')]);

figure;
plot(fwhm(:,:),'LineWidth',2); hold on; plot(durations(:),'LineWidth',2);
xlabel('repetition #','Color','k'); ylabel('Time from onset (sec)','Color','k'); legend({'Thr onset' 'Thr offset' 'durations'});
set(gca,'FontSize',32);
xlim([0 numel(durations)+1]);

figure;
ph = plot(durations(:),fwhm(:,:),'o','MarkerSize',14);
ph(1).MarkerFaceColor = 'b';
ph(2).MarkerFaceColor = 'r';
set(gca,'FontSize',32);
axis tight;

figure;%diff
ph = plot(durations(:),diff(fwhm(:,:)')','o','MarkerSize',14,'MarkerFaceColor','k','MarkerEdgeColor','none');
xlabel('Phrase duration (sec)','Color','k'); ylabel('Ca^{+2} width (sec)','Color','k');
set(gca,'FontSize',32);
axis tight;

figure;
plot(results(Day).Max(cellnum,:),'ro','MarkerSize',14,'MarkerFaceColor','r'); set(gca,'XTick',1:41); set(gca,'XTickLabel',syllables);
xlabel('syllables'); ylabel('0.9 quantile df/f');
xtickangle(90);
set(gca,'FontSize',18);


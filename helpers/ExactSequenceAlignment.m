function elements = ExactSequenceAlignment(Day,sylnum)
%%
% Day = '2017_06_28';
% sylnum = 8;
tmpthr=0;
startwith = 1;
win_edge_size = 0.05;
Nsyls_in_win = 10;
min_gap = 0.001;
min_syl = 0.005;
%% Folders that contain data
% Folders on laptop:
bird_name = 'lrb85315';
bird_folder_name = 'lrb853_15';
template_file = 'lrb85315template';
annotation_file = 'lrb85315auto_annotation5_fix'; %'lrb85315auto_annotation5_fix';
CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

laptop_manualROI_analyses_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name '/ManualROIs/PhraseLockedSpikeTiming'];
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));

%%


cd (laptop_manualROI_folder);
load(template_file,'templates');
syllables = [[templates.wavs.segType] -1 102 103];
n_syllables = numel(syllables);
freq_min = 300; freq_max = 8000;
colors = distinguishable_colors(n_syllables,'w');
load(annotation_file,'keys','elements');
ord = [];
dates = [];
issyl = [];
phrase_durations = [];
for i = 1:numel(keys)
    tokens = regexp(keys{i},'_','split');
    ord = [ord; str2num(tokens{2})];
    dates = [dates; char(join(tokens(3:5),'_'))];
    issyl = [issyl; ismember(sylnum,elements{i}.segType)];
    if ismember(sylnum,elements{i}.segType)
        phrases = return_phrase_times(elements{i});
        phrase_durations = [phrase_durations; ...
            max(phrases.phraseFileEndTimes(phrases.phraseType == sylnum)-phrases.phraseFileStartTimes(phrases.phraseType == sylnum))];
    else
        phrase_durations = [phrase_durations;0];
    end
end
[locs,indx] = sort(ord);
elements = elements(indx);
keys = keys(indx);
dates = dates(indx,:);
unique_dates = datestr(setdiff(unique(datenum(dates)),[736804]),'yyyy_mm_dd'); %does not include 04/19-21th (remove for other birds)
%%

res = [];
fmax = 8000;
files_idx = find(datenum(dates) == datenum(Day));
FILES = keys(files_idx);
maxwidth = ceil(max(phrase_durations(files_idx))*1000);
all_syl_durations = [];
cell_syl_durations = {};
hits = [];
for file_cnt = startwith:numel(files_idx)
    filename = keys{files_idx(file_cnt)};
    [y,fs] = audioread(fullfile(laptop_wav_folder,filename));
    [S,F,T,P] = spectrogram((y/(sqrt(mean(y.^2)))),220,220-44,512,fs);%,'reassigned');
    phrases = return_phrase_times(elements{files_idx(file_cnt)});
    locs = find(phrases.phraseType == sylnum);
    for cnt = 1:numel(locs)
        tonset = phrases.phraseFileStartTimes(locs(cnt));
        toffset = phrases.phraseFileEndTimes(locs(cnt));
        syl_idx = find(elements{files_idx(file_cnt)}.segFileStartTimes >= tonset & ...
            elements{files_idx(file_cnt)}.segFileStartTimes < toffset);
        %display(elements{files_idx(file_cnt)}.segType(syl_idx(1)));
        display([num2str(numel(syl_idx)) ' syllable']);
        start_syl = 1;
        while (start_syl <= numel(syl_idx))
            t1 = elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(start_syl))-win_edge_size;
            t2 = elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(min(start_syl+Nsyls_in_win-1,numel(syl_idx))))+win_edge_size;
            h_temp = figure('Position',[72         392        1138         313]); 
            ax_temp = axes;
            plot(T(T >= t1 & T<=t2),log(sum(abs(S(F<fmax & F>0,T >= t1 & T<=t2)))));
            xlim([t1 t2]);
            title(['File #' num2str(file_cnt) ' of ' num2str(numel(files_idx))]);
            if (tmpthr == 0)
                tmpthr = quantile(log(sum(abs(S(F<fmax & F>0,T >= t1 & T<=t2)))),0.1);
            end
            h_line = imline(ax_temp,[t1 tmpthr; t2 tmpthr]);
            pause;
            thr = (ax_temp.Children(1).Children(1).YData + ax_temp.Children(1).Children(2).YData)/2;
            tmpthr = thr;
            %hgclose(h_temp);
            [on_times,off_times] = syllable_envelope(log(sum(abs(S(F<fmax & F>0,T >= t1 & T<=t2)))),T(T >= t1 & T<=t2),thr,min_gap,min_syl);
            % prepare mock syllables for auto positioning
            mock_offs = off_times(off_times > on_times(1));
            mock_ons = on_times(on_times < mock_offs(end));
            mock_centers = (mock_offs+mock_ons)/2;
            %
            
            hf=figure('Position',[72         26        1138         313]); 
            ax = axes;
            %set(ax,'ButtonDownFcn','hselect_Callback');
        %hpanel=uipanel('position',[0 .05 2 .95]);
        %hscrollbar=uicontrol('style','slider','units','normalized','position',[0 0 1 .05],'callback',@hscroll_Callback);
        %axes('parent',hpanel,'outerposition',[.25 0 .5 1])
            imagesc(ax,T(T >= t1 & T<=t2),F(F<fmax),abs(S(F<fmax,T >= t1 & T<=t2))); colormap(1-gray); caxis([0 15]); xlim([t1 t2]);
            xlim([t1 t2]);
            set(gca,'YDir','normal');
            hold on; 
            for line_cnt = 1:numel(on_times)
                line([on_times(line_cnt) on_times(line_cnt)],[0 fmax],'Color',[0 0.7 0]);
            end
            for line_cnt = 1:numel(off_times)
                line([off_times(line_cnt) off_times(line_cnt)],[0 fmax],'Color',[0.7 0 0]);
            end
            hs=[];
            h=0;
            for syl_num = start_syl:min(start_syl+Nsyls_in_win-1,numel(syl_idx))
                h=imrect(ax,[elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) ...
                     0 elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num)) - ...
                     elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) 8000]);
                 hs=[hs; h];
            end
            set(hf,'WindowbuttonDownFcn',@clickcallback)
            set(hf,'KeyPressFcn',@(h_obj,evt) keystroke(h_obj,evt));
            drawnow;
            flag = 0;
            while flag == 0
                pause(0.5)
            end
           
            
            for syl_num = start_syl:min(start_syl+numel(hs)-1,numel(syl_idx))
                syl_cnt = syl_num - start_syl+1;
                h = get(hs(syl_cnt),'Children');
                maxx = h(1).XData;
                minx = h(3).XData;
                elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) = minx;
                elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num)) = maxx;
            end
            '*';
            hgclose(hf);
            hgclose(h_temp);
            start_syl = min(start_syl+numel(hs)-1,numel(syl_idx))+1;
        end
       % row = zeros(1,maxwidth);
         syl_durations = elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx) - ...
             elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx);
%         cell_syl_durations = {cell_syl_durations{:} syl_durations};
%         all_syl_durations = [all_syl_durations syl_durations];
%         for syl_num = 1:numel(syl_idx)
%             bins = max(floor(1000*(elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num))-tonset))+1,1):min(ceil(1000*(elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num))-tonset))-1,maxwidth);
%             row(bins) = 1;
%         end
%         res = [res; row];
        hits = [hits; files_idx(file_cnt) locs(cnt) numel(syl_idx) min(syl_durations) max(syl_durations)];
    end
end

    function clickcallback(obj,evt)
        hd = get(obj,'Children');
        persistent chk
        if isempty(chk)
              chk = 1;
              pause(0.5); %Add a delay to distinguish single click from a double click
              if chk == 1
                  %fprintf(1,'\nI am doing a single-click.\n\n');
                  chk = [];
              end
        else
              chk = [];
              mouse_loc = get(hd(end),'CurrentPoint');
              xpos = mouse_loc(1,1); ypos = mouse_loc(1,2);
              if (xpos > t1 & xpos < t2)
                  hitsyl = [];
                  for syl_num = start_syl:min(start_syl+numel(hs)-1,numel(syl_idx))
                        syl_cnt = syl_num - start_syl+1;
                        h = get(hs(syl_cnt),'Children');
                        maxx = h(1).XData;
                        minx = h(3).XData;
                        if (xpos > minx & xpos < maxx)
                            hitsyl = sylnum;
                            rectpos = getPosition(hs(syl_cnt));
                            rec_center = rectpos(1) + rectpos(3)/2;
                            mock_dist = abs(mock_centers - rec_center);
                            mock_loc = min(find(mock_dist == min(mock_dist)));
                            setPosition(hs(syl_cnt),[mock_ons(mock_loc) rectpos(2) (mock_offs(mock_loc)-mock_ons(mock_loc)) rectpos(4)]);
                        end
                  end
              end
              %fprintf(1,'\nI am doing a double-click.\n\n');
        end
    end

    function keystroke(h_obj,evt)
        
        disp(evt.Key);
        hd = get(h_obj,'Children');
        mouse_loc = get(hd(end),'CurrentPoint');
        xpos = mouse_loc(1,1); ypos = mouse_loc(1,2);
        disp([xpos ypos]);
        switch evt.Key
            case 'u'
                thr = (ax_temp.Children(1).Children(1).YData + ax_temp.Children(1).Children(2).YData)/2;
                tmpthr = thr;
                for syl_num = start_syl:min(start_syl+numel(hs)-1,numel(syl_idx))
                    syl_cnt = syl_num - start_syl+1;
                    h = get(hs(syl_cnt),'Children');
                    maxx = h(1).XData;
                    minx = h(3).XData;
                    elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) = minx;
                    elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num)) = maxx;
                end
                [on_times,off_times] = syllable_envelope(log(sum(abs(S(F<fmax & F>0,T >= t1 & T<=t2)))),T(T >= t1 & T<=t2),thr,min_gap,min_syl);
                % prepare mock syllables for auto positioning
                mock_offs = off_times(off_times > on_times(1));
                mock_ons = on_times(on_times < mock_offs(end));
                mock_centers = (mock_offs+mock_ons)/2;
                axes(ax); hold off;
                imagesc(ax,T(T >= t1 & T<=t2),F(F<fmax),abs(S(F<fmax,T >= t1 & T<=t2))); colormap(1-gray); caxis([0 15]); xlim([t1 t2]);
                xlim([t1 t2]);
                set(gca,'YDir','normal');
                hold on; 
                for line_cnt = 1:numel(on_times)
                    line([on_times(line_cnt) on_times(line_cnt)],[0 fmax],'Color',[0 0.7 0]);
                end
                for line_cnt = 1:numel(off_times)
                    line([off_times(line_cnt) off_times(line_cnt)],[0 fmax],'Color',[0.7 0 0]);
                end
                hs=[];
                h=0;
                for syl_num = start_syl:min(start_syl+Nsyls_in_win-1,numel(syl_idx))
                    h=imrect(ax,[elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) ...
                         0 elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num)) - ...
                         elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) 8000]);
                     hs=[hs; h];
                end
                set(hf,'WindowbuttonDownFcn',@clickcallback)
                set(hf,'KeyPressFcn',@(h_obj,evt) keystroke(h_obj,evt));
                drawnow;
                
            case 'a'
                for syl_num = start_syl:min(start_syl+numel(hs)-1,numel(syl_idx))
                    syl_cnt = syl_num - start_syl+1;
                    h = get(hs(syl_cnt),'Children');
                    maxx = h(1).XData;
                    minx = h(3).XData;
                    elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) = minx;
                    elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num)) = maxx;
                end
%                 rectpos = getPosition(hs(1)); minx = rectpos(1);
%                 rectpos = getPosition(hs(end)); maxx = rectpos(1)+rectpos(3);
                rect = getrect(ax);
                new_rec_minx = rect(1); new_rec_maxx = rect(1)+rect(3);
%                 if (new_rec_minx < maxx & new_rec_minx > minx & new_rec_maxx < maxx & new_rec_maxx > minx)
                hitsyl = 0;
                rectpos = getPosition(hs(1));
                if (new_rec_maxx < rectpos(1))
                    hitsyl = 1;
                    h=imrect(ax,[rect(1) 0 rect(3) 8000]);
                    hs = [h; hs];
                end
                rectpos = getPosition(hs(end));
                if (new_rec_minx > rectpos(1)+rectpos(3))
                    hitsyl = numel(hs)+1;
                    h=imrect(ax,[rect(1) 0 rect(3) 8000]);
                    hs = [hs; h];
                end
                if (hitsyl == 0)
                    for syl_num = start_syl:min(start_syl+numel(hs)-2,numel(syl_idx))
                        syl_cnt = syl_num - start_syl+1;
                        rectpos = getPosition(hs(syl_cnt));
                        minx = rectpos(1)+rectpos(3);
                        rectpos = getPosition(hs(syl_cnt+1));
                        maxx = rectpos(1);
                        if (new_rec_minx > minx & new_rec_maxx < maxx)
                            hitsyl = syl_cnt;
                            h=imrect(ax,[rect(1) 0 rect(3) 8000]);
                            hs = [hs(1:hitsyl); h; hs(hitsyl+1:end)];
                        end
                    end
                end
                
                if (hitsyl ~= 0)
                    syl_before = max(find(elements{files_idx(file_cnt)}.segFileStartTimes < rect(1)));
                    elements{files_idx(file_cnt)}.segFileStartTimes = ...
                        [elements{files_idx(file_cnt)}.segFileStartTimes(1:syl_before) rect(1) ...
                        elements{files_idx(file_cnt)}.segFileStartTimes(syl_before+1:end)];
                    elements{files_idx(file_cnt)}.segAbsStartTimes = ...
                        [elements{files_idx(file_cnt)}.segAbsStartTimes(1:syl_before) rect(1) + ...
                        elements{files_idx(file_cnt)}.segAbsStartTimes(syl_before) - elements{files_idx(file_cnt)}.segFileStartTimes(syl_before) ...
                        elements{files_idx(file_cnt)}.segAbsStartTimes(syl_before+1:end)];
                    elements{files_idx(file_cnt)}.segFileEndTimes = ...
                        [elements{files_idx(file_cnt)}.segFileEndTimes(1:syl_before) rect(1)+rect(3) ...
                        elements{files_idx(file_cnt)}.segFileEndTimes(syl_before+1:end)];
                    elements{files_idx(file_cnt)}.segType = ...
                        [elements{files_idx(file_cnt)}.segType(1:syl_before); sylnum; ...
                        elements{files_idx(file_cnt)}.segType(syl_before+1:end)];
                    phrases = return_phrase_times(elements{files_idx(file_cnt)});
                    tonset = phrases.phraseFileStartTimes(locs(cnt));
                    toffset = phrases.phraseFileEndTimes(locs(cnt));
                    syl_idx = find(elements{files_idx(file_cnt)}.segFileStartTimes >= tonset & ...
                        elements{files_idx(file_cnt)}.segFileStartTimes < toffset);
                    
                end
%         end
            case 'n'
                flag = 1;
            case 'd'
                if (xpos > t1 & xpos < t2)
                  hitsyl = [];
                  for syl_num = start_syl:min(start_syl+numel(hs)-1,numel(syl_idx))
                        syl_cnt = syl_num - start_syl+1;
                        rectpos = getPosition(hs(syl_cnt));
                        maxx = rectpos(1)+rectpos(3);
                        minx = rectpos(1);
                        if (xpos > minx & xpos < maxx)
                            hitsyl = sylnum;
                            elements{files_idx(file_cnt)}.segFileStartTimes(syl_idx(syl_num)) = [];
                            elements{files_idx(file_cnt)}.segFileEndTimes(syl_idx(syl_num)) = [];
                            elements{files_idx(file_cnt)}.segAbsStartTimes(syl_idx(syl_num)) = [];
                            elements{files_idx(file_cnt)}.segType(syl_idx(syl_num)) = [];
                            phrases = return_phrase_times(elements{files_idx(file_cnt)});
                            tonset = phrases.phraseFileStartTimes(locs(cnt));
                            toffset = phrases.phraseFileEndTimes(locs(cnt));
                            syl_idx = find(elements{files_idx(file_cnt)}.segFileStartTimes >= tonset & ...
                                elements{files_idx(file_cnt)}.segFileStartTimes < toffset);
                            
                            delete(hs(syl_cnt));
                            hs(syl_cnt) = [];
                            break;
                        end
                  end
              end
                
        end
    end
end
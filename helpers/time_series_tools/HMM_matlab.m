%%
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/pmtk3'),'-end');
nstates = 2;
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata_'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_'};
bird_params = bird2_params;

bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
file_prefix = bird_params{5}; 

CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

%%

cd (laptop_manualROI_folder);

load(annotation_file);
ord = [];
dates = [];
for i = 1:numel(keys)
    tokens = regexp(keys{i},'_','split');
    ord = [ord; str2num(tokens{2})];
    dates = [dates; char(join(tokens(3:5),'_'))];
end
[locs,indx] = sort(ord);
elements = elements(indx);
keys = keys(indx);
dates = dates(indx,:);
unique_dates = datestr(unique(datenum(dates)),'yyyy_mm_dd'); %datestr(setdiff(unique(datenum(dates)),736804),'yyyy_mm_dd'); %does not include 04/19th (remove for other birds)
%%
for Day_num = 1: size(unique_dates,1)
    Day = unique_dates(Day_num,:);
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    FILES = dir([file_prefix bird_name '*.mat']);
    %FILES = dir('NonoverlapBaseROIdata*.mat');
    FILES = {FILES.name};
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        load(fname);
        map_states = zeros(size(dff));
        for roi_n = 1:size(dff,1)
            sig = dff(roi_n,:);
            clear sigma;
            try
                mu = [quantile(sig,0.05) max(sig)];
                sigma(1,1,:) = [0.05*ones(1,nstates)];
                CPD = condGaussCpdCreate(mu, sigma);
                [model, loglikHist] = hmmFit(sig, nstates, 'gauss','pi0',[zeros(1,nstates-1) 1],'emission0',CPD);%,'transPrior',[10 1;10 3],'piPrior',[100 1]);                    
                map_path = hmmMap(model, sig)-1;
                map_path = abs(median(map_path)-map_path);
                map_states(roi_n,:) = map_path;
            catch em
               '-'; 
               display([Day_num fnum roi_n]);
            end
        end
        outname = ['GaussHmmMAP_' fname];
        save(outname,'map_states');
    end
end
%% group-level RSA effect
clear;clc
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_statisticsScripts\circstat-matlab-master\circstat-matlab-master
dataDir = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_noneEEG_preProcessedBaselines\Cardiac_Breathing';

% Get a list of all folders in the directory
contents = dir(dataDir);
folders = contents([contents.isdir]);
subj_folders = folders(~ismember({folders.name}, {'.', '..','03','06','08','20','24','VTD_testRecordings'})); %03,06,20,24 are excluded from breathing analyses

% Initialize an empty cell array to store the numeric parts as strings
number_strings = {};

% Loop through the folder list to extract the numeric parts
for i = 1:length(subj_folders)
    folder_name = subj_folders(i).name; % Get the folder name
    % Use regular expression to extract the numeric part
    num_str = regexp(folder_name, '\d+', 'match');
    if ~isempty(num_str)
        % Store the numeric part as a string with leading zeros (if necessary)
        formatted_num_str = sprintf('%02d', str2double(num_str{1}));
        number_strings = [number_strings, formatted_num_str];
    end
end

for subj=1:size(number_strings,2)
    subjData = [dataDir '\' number_strings{1,subj}];
    rsaData(subj) = load([subjData '\rsa_effectBaseline_fig2B.mat']);
end

% cat1_gMean = mean([rsaData.cat1_mean]); cat1_gSEM = std([rsaData.cat1_mean])/sqrt(length([rsaData.cat1_mean]));
% cat2_gMean = mean([rsaData.cat2_mean]); cat2_gSEM = std([rsaData.cat2_mean])/sqrt(length([rsaData.cat2_mean]));
% cat3_gMean = mean([rsaData.cat3_mean]); cat3_gSEM = std([rsaData.cat3_mean])/sqrt(length([rsaData.cat3_mean]));
% cat4_gMean = mean([rsaData.cat4_mean]); cat4_gSEM = std([rsaData.cat4_mean])/sqrt(length([rsaData.cat4_mean]));
% cat5_gMean = mean([rsaData.cat5_mean]); cat5_gSEM = std([rsaData.cat5_mean])/sqrt(length([rsaData.cat5_mean]));
% cat6_gMean = mean([rsaData.cat6_mean]); cat6_gSEM = std([rsaData.cat6_mean])/sqrt(length([rsaData.cat6_mean]));
% cat7_gMean = mean([rsaData.cat7_mean]); cat7_gSEM = std([rsaData.cat7_mean])/sqrt(length([rsaData.cat7_mean]));
% cat8_gMean = mean([rsaData.cat8_mean]); cat8_gSEM = std([rsaData.cat8_mean])/sqrt(length([rsaData.cat8_mean]));

%% Bar graph representation
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\ViolinPlot\daboxplot
figure;
hold on;
% mean_hr = [[rsaData.cat1_mean] [rsaData.cat2_mean] [rsaData.cat3_mean] [rsaData.cat4_mean] [rsaData.cat5_mean] [rsaData.cat6_mean] [rsaData.cat7_mean] [rsaData.cat8_mean]];
% se_hr = [0, std(mean4590)/sqrt(length(mean4590)), std(mean90135)/sqrt(length(mean90135)), std(mean135180)/sqrt(length(mean135180)), std(mean180225)/sqrt(length(mean180225)), std(mean225270)/sqrt(length(mean225270)), std(mean270315)/sqrt(length(mean270315)), std(mean315360)/sqrt(length(mean315360))];

% meanzeros = zeros(23,1);
colors = repmat([0.5 0.5 0.5],8,1);
daboxplot([[rsaData.cat1_mean]' [rsaData.cat2_mean]' [rsaData.cat3_mean]' [rsaData.cat4_mean]' [rsaData.cat5_mean]' [rsaData.cat6_mean]' [rsaData.cat7_mean]' [rsaData.cat8_mean]'],'colors',colors,'outliers',1);
set(gca,'FontSize',36);
ylabel('Delta IBI');
xlabel('Breathing Phase');
xticklabels({'0-45','45-90','90-135','135-180','180-225','225-270','270-315','315-360'})
title('Changes in IBI across the respiratory cycle')
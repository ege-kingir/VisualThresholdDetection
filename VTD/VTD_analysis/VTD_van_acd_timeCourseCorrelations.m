%% looking into the between-subject correlations between the VAN time courses and the cardiac deceleration parameters
clear;clc

%%% CARDIAC DECEL PARAMS (for Hit-Miss):
%   a) ACD (total cardiac deceleration): (mean_acd_s0_hit - mean_acd_s2_hit) - (mean_acd_s0_miss - mean_acd_s2_miss);
%   b) delta B: ([rsa_ibi_modelData_BBB.B_hit_s0]'-[rsa_ibi_modelData_BBB.B_hit_s2]')-([rsa_ibi_modelData_BBB.B_miss_s0]'-[rsa_ibi_modelData_BBB.B_miss_s2]');
%   c) A: [rsa_ibi_modelData.A_hit]' - [rsa_ibi_modelData.A_miss]'
%% load the downsampled Hit vs. Miss stimulus contrast
stimConts = load("S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing_newFilters\hitAndMiss_stimContrasts.mat");
hitConts = stimConts.hit_stimContrast;
hit_dsConts = stimConts.hit_ds_stimContrast;
missConts = stimConts.miss_stimContrast;
avgConts = mean([hit_dsConts missConts],2);
contDiff = hit_dsConts - missConts;
%% load the van time course data ([nSubj x nTimePoints])
% clear;clc
CNV_data = load("S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing_newFilters\van_avg_timeCourses_hitsDS.mat");
CNV_dataAll = CNV_data.van_avgTimeCourses;
CNV_dataHit = CNV_data.van_avgTimeCoursesHit;
CNV_dataMiss = CNV_data.van_avgTimeCoursesMiss;
timeCNV = -100:2:1000;

%% specify for VAN (should be the central electrodes of Cz and CPz -- because it is the most prominent there)
CNV_VAN = CNV_dataAll(:,find(timeCNV==-100):find(timeCNV==1000));
CNV_baseline = CNV_dataAll(:,find(timeCNV==-100):find(timeCNV==0));
CNV_VAN = CNV_VAN-mean(CNV_baseline,2);
CNV_z = zscore(CNV_VAN, 0, 1);           % z-score across subjects for each time point

%% specify for VAN -- only Hits
CNV_VANhit = CNV_dataHit(:,find(timeCNV==-100):find(timeCNV==1000));
CNV_baselineHit = CNV_dataHit(:,find(timeCNV==-100):find(timeCNV==0));
CNV_VANhit = CNV_VANhit - mean(CNV_baselineHit,2);
CNV_zHit = zscore(CNV_VANhit,0,1);

%% specify for VAN -- only Misses
CNV_VANmiss = CNV_dataMiss(:,find(timeCNV==-100):find(timeCNV==1000));
CNV_baselineMiss = CNV_dataMiss(:,find(timeCNV==-100):find(timeCNV==0));
CNV_VANmiss = CNV_VANmiss - mean(CNV_baselineMiss,2);
CNV_zMiss = zscore(CNV_VANmiss,0,1);

CNV_VANdiff = CNV_VANhit - CNV_VANmiss;

%% load the RSA model parameter data
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_statisticsScripts\circstat-matlab-master\circstat-matlab-master
dataDir = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData';

% Get a list of all folders in the directory
contents = dir(dataDir);
folders = contents([contents.isdir]);
subj_folders = folders(~ismember({folders.name}, {'.', '..','VTD03','VTD06','VTD08','VTD20','VTD24','VTD_testRecordings'}));

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
    subjData = [dataDir '\VTD' number_strings{1,subj}];
    rsa_ibi_modelData(subj) = load([subjData '\rsa_ibi_sineModel.mat']);
end

for subj=1:size(number_strings,2) % BBB=beat by beat
    subjData = [dataDir '\VTD' number_strings{1,subj}];
    rsa_ibi_modelData_BBB(subj) = load([subjData '\rsa_ibi_sineModel_beatByBeat.mat']);
end

%% load the ACD data (from ECG signals)
for subj=1:size(number_strings,2)
    subjData = [dataDir '\VTD' number_strings{1,subj}];
    acd(subj) = load([subjData '\acd_preStim_brPhase.mat']);
end
nSubj=23;
for s=1:nSubj
    mean_acd_s2_hit(s,1) = mean(acd(s).acd_s2_hit,'omitnan'); mean_acd_s2_miss(s,1) = mean(acd(s).acd_s2_miss,'omitnan');
    mean_acd_s1_hit(s,1) = mean(acd(s).acd_s1_hit,'omitnan'); mean_acd_s1_miss(s,1) = mean(acd(s).acd_s1_miss,'omitnan');
    mean_acd_s0_hit(s,1) = mean(acd(s).acd_s0_hit,'omitnan'); mean_acd_s0_miss(s,1) = mean(acd(s).acd_s0_miss,'omitnan');
%     mean_acd_sp1_hit(s,1) = mean(acd(s).acd_s_p1_hit,'omitnan'); mean_acd_sp1_miss(s,1) = mean(acd(s).acd_s_p1_miss,'omitnan');
    
    sem_acd_s2_hit(s,1) = std(acd(s).acd_s2_hit,'omitnan')/sqrt(size(acd(s).acd_s2_hit,2)); sem_acd_s2_miss(s,1) = std(acd(s).acd_s2_miss,'omitnan')/sqrt(size(acd(s).acd_s2_miss,2));
    sem_acd_s1_hit(s,1) = std(acd(s).acd_s1_hit,'omitnan')/sqrt(size(acd(s).acd_s1_hit,2)); sem_acd_s1_miss(s,1) = std(acd(s).acd_s1_miss,'omitnan')/sqrt(size(acd(s).acd_s1_miss,2));
    sem_acd_s0_hit(s,1) = std(acd(s).acd_s0_hit,'omitnan')/sqrt(size(acd(s).acd_s0_hit,2)); sem_acd_s0_miss(s,1) = std(acd(s).acd_s0_miss,'omitnan')/sqrt(size(acd(s).acd_s0_miss,2));
%     sem_acd_sp1_hit(s,1) = std(acd(s).acd_s_p1_hit,'omitnan')/sqrt(size(acd(s).acd_s_p1_hit,2)); sem_acd_sp1_miss(s,1) = std(acd(s).acd_s_p1_miss,'omitnan')/sqrt(size(acd(s).acd_s_p1_miss,2));
    mean_acd_s2(s,1) = mean(acd(s).acd_s2,'omitnan');
    mean_acd_s1(s,1) = mean(acd(s).acd_s1,'omitnan');
    mean_acd_s0(s,1) = mean(acd(s).acd_s0,'omitnan');

    sem_acd_s2(s,1) = std(acd(s).acd_s2,'omitnan')/sqrt(size(acd(s).acd_s2,2));
    sem_acd_s1(s,1) = std(acd(s).acd_s1,'omitnan')/sqrt(size(acd(s).acd_s1,2));
    sem_acd_s0(s,1) = std(acd(s).acd_s0,'omitnan')/sqrt(size(acd(s).acd_s0,2));
end

%% Input
deltaOffset = mean_acd_s0-mean_acd_s2;

%% get the correlation between the significant VAN window (significant difference between Hit and Miss = 170-318 ms // dsHits = 188-318 ms)
vanSign = CNV_VAN(:,find(timeCNV==188):find(timeCNV==318));
x_c = mean(vanSign,2);
y_c = deltaOffset;
% Fit a linear regression model (degree 1)
p = polyfit(x_c, y_c, 1);

% p(1) is the slope, p(2) is the intercept
slope = p(1);
intercept = p(2);

% Display the results
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);

figure;
% To plot the data and the fitted line
y_fit = polyval(p, x_c);
plot(x_c, y_c, 'o','MarkerSize',20,'Color','k','LineWidth',5);         % Plot original data points
hold on;
plot(x_c, y_fit, '-','LineWidth',10,'Color',[0.5 0 0]);     % Plot the fitted line
hold off;
xlabel('VAN');
ylabel('ACD');
title('VAN vs ACD: All');
legend('Data', 'Fitted Line');

% Add text to the northwest corner
text_loc_x = min(x_c); % Use the minimum x value for positioning
text_loc_y = max(y_c); % Use the maximum y value for positioning

% Calculate correlation coefficient and p-value
mdl_huber = fitlm(x_c, y_c, 'RobustOpts', 'huber');
rAll = mdl_huber.Rsquared.Ordinary;
rAll = sqrt(rAll)

[r p] = corr(x_c, y_c,'type','Spearman')

set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
text(text_loc_x, text_loc_y, ...
    sprintf('r = %.2f\np = %.4f', mdl_huber.Coefficients.Estimate(2), mdl_huber.Coefficients.pValue(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline', ...
    'FontSize', 48, 'BackgroundColor', 'white', 'EdgeColor', 'black');
set(gca, 'FontSize', 36);

%% correlation analysis based on only Hit trials
% Inputs
vanHit = CNV_VANhit(:,find(timeCNV==188):find(timeCNV==318));
deltaOffsetH = [rsa_ibi_modelData.A_hit]'; %mean_acd_s0_hit-mean_acd_s2_hit %[rsa_ibi_modelData_BBB.B_hit_s0]'-[rsa_ibi_modelData_BBB.B_hit_s2]'
x_c = mean(vanHit,2);
y_c = deltaOffsetH;
% Fit a linear regression model (degree 1)
p = polyfit(x_c, y_c, 1);

% p(1) is the slope, p(2) is the intercept
slope = p(1);
intercept = p(2);

% Display the results
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);

figure;
% To plot the data and the fitted line
y_fit = polyval(p, x_c);
plot(x_c, y_c, 'o','MarkerSize',20,'Color','k','LineWidth',5);         % Plot original data points
hold on;
plot(x_c, y_fit, '-','LineWidth',10,'Color',[0.5 0 0]);     % Plot the fitted line
hold off;
xlabel('VAN');
ylabel('ACD');
title('VAN vs ACD: Only Hits');
legend('Data', 'Fitted Line');

% Add text to the northwest corner
text_loc_x = min(x_c); % Use the minimum x value for positioning
text_loc_y = max(y_c); % Use the maximum y value for positioning

% Calculate correlation coefficient and p-value
mdl_huber = fitlm(x_c, y_c, 'RobustOpts', 'huber')
rHit = mdl_huber.Rsquared.Ordinary;
rHit = sqrt(rHit)

% get the same rho from corr()
[rHit pHit] = corr(x_c, y_c,'type','Spearman')

set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
text(text_loc_x, text_loc_y, ...
    sprintf('r = %.2f\np = %.4f', mdl_huber.Coefficients.Estimate(2), mdl_huber.Coefficients.pValue(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline', ...
    'FontSize', 48, 'BackgroundColor', 'white', 'EdgeColor', 'black');
set(gca, 'FontSize', 36);

%% is the across-subject correlation of VAN~ACD (or delta B), time-window specific???
deltaOffset = mean_acd_s0_hit - mean_acd_s2_hit;
cardiac_z = zscore(deltaOffset);

nSubjects = size(CNV_VANhit, 1);
nTimePoints = size(CNV_VANhit, 2);
nPermutations = 5000;
cluster_alpha = 0.01;   % threshold for forming clusters
perm_alpha = 0.01;      % significance level after correction

% Step 1: Compute observed Spearman correlations
r_obs = zeros(1, nTimePoints);
p_obs = zeros(1, nTimePoints);
for t = 1:nTimePoints
    [r_obs(t), p_obs(t)] = corr(CNV_VANhit(:,t), deltaOffset, 'Type', 'Spearman');
end

% Step 2: Identify clusters of significant correlations
cluster_inds = bwconncomp(p_obs < cluster_alpha);

% Sum of absolute r-values within each cluster
cluster_r_obs = zeros(1, cluster_inds.NumObjects);
for i = 1:cluster_inds.NumObjects
    cluster_r_obs(i) = sum(abs(r_obs(cluster_inds.PixelIdxList{i})));
end

% Step 3: Permutation testing
max_cluster_r_perm = zeros(1, nPermutations);
for p = 1:nPermutations
    permuted_card = deltaOffset(randperm(nSubjects));
    r_perm = zeros(1, nTimePoints);
    for t = 1:nTimePoints
        r_perm(t) = corr(CNV_VANhit(:,t), permuted_card, 'Type', 'Spearman');
    end
    % Find clusters in permuted data
    cluster_perm = bwconncomp(abs(r_perm) > quantile(abs(r_obs), 1 - cluster_alpha));
    if cluster_perm.NumObjects > 0
        cluster_sums = zeros(1, cluster_perm.NumObjects);
        for i = 1:cluster_perm.NumObjects
            cluster_sums(i) = sum(abs(r_perm(cluster_perm.PixelIdxList{i})));
        end
        max_cluster_r_perm(p) = max(cluster_sums);
    end
end

% Step 4: Compare observed clusters to permutation distribution
significant_clusters = [];
for i = 1:length(cluster_r_obs)
    if cluster_r_obs(i) > prctile(max_cluster_r_perm, 100 * (1 - perm_alpha))
        significant_clusters = [significant_clusters i];
    end
end

% % Step 5: Plot
time = linspace(-0.1,1, nTimePoints);  % adjust if needed
figure; hold on;
plot(time, r_obs, 'k','LineWidth',5);
yline(0, '--', 'Color', [0.5 0.5 0.5]);

for i = significant_clusters
    idx = cluster_inds.PixelIdxList{i};
    plot(time(idx), r_obs(idx), 'b', 'LineWidth', 10);
end

% allSignTime = time(idx);

xlabel('Time (s)');
ylabel('Spearman r');
title('Occipital electrodes: VAN ~ ACD');
legend('Observed r', '', 'Significant cluster');
set(gca,'FontSize',36);

%% correlation analysis based on only Miss trials
% Inputs
vanMiss = CNV_VANmiss(:,find(timeCNV==188):find(timeCNV==318));
deltaOffsetM = [rsa_ibi_modelData.A_miss]'; %[rsa_ibi_modelData_BBB.B_miss_s0]'-[rsa_ibi_modelData_BBB.B_miss_s2]'; %mean_acd_s0_miss-mean_acd_s2_miss
x_c = mean(vanMiss,2);
y_c = deltaOffsetM;
% Fit a linear regression model (degree 1)
p = polyfit(x_c, y_c, 1);

% p(1) is the slope, p(2) is the intercept
slope = p(1);
intercept = p(2);

% Display the results
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);

figure;
% To plot the data and the fitted line
y_fit = polyval(p, x_c);
plot(x_c, y_c, 'o','MarkerSize',20,'Color','k','LineWidth',5);         % Plot original data points
hold on;
plot(x_c, y_fit, '-','LineWidth',10,'Color',[0.5 0 0]);     % Plot the fitted line
hold off;
xlabel('VAN');
ylabel('ACD');
title('VAN vs ACD: Only Misses');
legend('Data', 'Fitted Line');

% Add text to the northwest corner
text_loc_x = min(x_c); % Use the minimum x value for positioning
text_loc_y = max(y_c); % Use the maximum y value for positioning

% Calculate correlation coefficient and p-value
mdl_huber = fitlm(x_c, y_c, 'RobustOpts', 'huber')
rMiss = mdl_huber.Rsquared.Ordinary;
rMiss = sqrt(rMiss)

% get the same rho from corr()
[rMiss pMiss] = corr(x_c, y_c,'type','Spearman')

set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
text(text_loc_x, text_loc_y, ...
    sprintf('r = %.2f\np = %.4f', mdl_huber.Coefficients.Estimate(2), mdl_huber.Coefficients.pValue(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline', ...
    'FontSize', 48, 'BackgroundColor', 'white', 'EdgeColor', 'black');
set(gca, 'FontSize', 36);

%% look at the statistical difference between the rHit and rMiss
jh = fitlm(mean(vanHit,2),mean(vanMiss,2),'RobustOpts','huber'); rJH = sqrt(jh.Rsquared.Ordinary);
jm = fitlm(mean(vanHit,2),deltaOffsetM,'RobustOpts','huber'); rJM = sqrt(jm.Rsquared.Ordinary);
kh = fitlm(mean(vanMiss,2),deltaOffsetH,'RobustOpts','huber'); rKH = sqrt(kh.Rsquared.Ordinary);
km = fitlm(deltaOffsetH,deltaOffsetM,'RobustOpts','huber'); rKM = sqrt(km.Rsquared.Ordinary);
n=23;

%% can the effect be explained by the stimulus contrast differences of each subject
vanSign = CNV_VANhit(:,find(timeCNV==188):find(timeCNV==318));
x_c = mean(vanSign,2);
y_c = hit_dsConts;
% Fit a linear regression model (degree 1)
p = polyfit(x_c, y_c, 1);

% p(1) is the slope, p(2) is the intercept
slope = p(1);
intercept = p(2);

% Display the results
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);

figure;
% To plot the data and the fitted line
y_fit = polyval(p, x_c);
plot(x_c, y_c, 'o','MarkerSize',20,'Color','k','LineWidth',5);         % Plot original data points
hold on;
plot(x_c, y_fit, '-','LineWidth',10,'Color',[0.5 0 0]);     % Plot the fitted line
hold off;
xlabel('VAN: Hit - Miss');
ylabel('Stim. Contrast: Hit - Miss');
title('VAN vs Contrast');
legend('Data', 'Fitted Line');

% Add text to the northwest corner
text_loc_x = min(x_c); % Use the minimum x value for positioning
text_loc_y = max(y_c); % Use the maximum y value for positioning

% Calculate correlation coefficient and p-value
% [r, p_value] = corr(x_c, y_c,'type','Spearman');
mdl_huber = fitlm(x_c, y_c, 'RobustOpts', 'huber')
rCont = mdl_huber.Rsquared.Ordinary;
rCont = sqrt(rCont)

set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
text(text_loc_x, text_loc_y, ...
    sprintf('slope = %.2f\np = %.4f', mdl_huber.Coefficients.Estimate(2), mdl_huber.Coefficients.pValue(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline', ...
    'FontSize', 48, 'BackgroundColor', 'white', 'EdgeColor', 'black');
set(gca, 'FontSize', 36);

%% Plot the VAN time courses from Hit and Miss trials
figure;
plot(time,mean(CNV_VANhit,1),'g','LineWidth',5);
hold on
plot(time,mean(CNV_VANmiss,1),'r','LineWidth',5);
hold on

x_start = allSignTime(1);
x_end = allSignTime(end);
y_limits = ylim;  % Use current y-axis limits

Add transparent gray rectangle
fill([x_start x_end x_end x_start], ...
     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
     [0.7 0.7 0.7], ...       % RGB color for gray
     'FaceAlpha', 0.3, ...    % Transparency
     'EdgeColor', 'none');    % No border

xlabel('Time (s)');
ylabel('VAN (microV)');
title('Occipital electrodes: VAN in Hit vs. Miss');
set(gca,'FontSize',36);
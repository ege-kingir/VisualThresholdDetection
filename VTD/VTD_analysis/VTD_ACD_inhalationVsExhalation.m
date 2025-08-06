%% %% VTD: Group-Level analysis of the cardiac deceleration params...
%%% Effect of time within anticipatory period
%%% Effect of Hit vs Miss
%%% Effect of Confidence (Low vs High)

%% Also includes linear correlation analyses
%%% Correlations between A and B params of RSA models
%%% Correlations between RSA estimation methods

%% Organize data
clear;clc
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_statisticsScripts\circstat-matlab-master\circstat-matlab-master
dataDir = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData';
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_rsaSinusoidalModel

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
    breathingData(subj) = load([subjData '\acd_inhalationVsExhalation.mat']); %if you want the approach where trials are categorized according to only stimulus onset, get "acd_inhalationVsExhalation_stimOnset"
end

for subj=1:size(number_strings,2)
    subjData = [dataDir '\VTD' number_strings{1,subj}];
    rsa_ibi_modelData(subj) = load([subjData '\rsa_ibi_sineModel.mat']);
end

for subj=1:size(number_strings,2) % BBB=beat by beat
    subjData = [dataDir '\VTD' number_strings{1,subj}];
    rsa_ibi_modelData_BBB(subj) = load([subjData '\rsa_ibi_sineModel_beatByBeat.mat']);
end

for subj=1:size(number_strings,2)
    subjData = [dataDir '\VTD' number_strings{1,subj}];
    acd(subj) = load([subjData '\acd_preStim_brPhase.mat']);
end
%% baseline RSA data
baseDir = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_noneEEG_preProcessedBaselines\Cardiac_Breathing';
for subj=1:size(number_strings,2)
    subjData = [baseDir '\' number_strings{1,subj}];
    baseRSA(subj) = load([subjData '\rsa_ibi_sineModelBaseline.mat']);
end

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\ViolinPlot\dabarplot
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\ViolinPlot\daboxplot

%% ACD progression plots -- Hit vs. Miss
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
var11 = mean_acd_s2_hit;
var12 = mean_acd_s1_hit;
var13 = mean_acd_s0_hit;
var21 = mean_acd_s2_miss;
var22 = mean_acd_s1_miss;
var23 = mean_acd_s0_miss;
var1 = mean([var11 var12 var13],2);
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];
var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_hit = std(var11)/sqrt(length(var11)); sem_minus2_miss = std(var21)/sqrt(length(var21));
sem_minus1_hit = std(var12)/sqrt(length(var12)); sem_minus1_miss = std(var22)/sqrt(length(var22));
sem_zero_hit = std(var13)/sqrt(length(var13)); sem_zero_miss = std(var23)/sqrt(length(var23));

%%% plot the IBI progression over time in Fast and Slow responses separately
figure
times = [1:3];

err1 = sem_minus2_hit;
err2 = sem_minus1_hit;
err3 = sem_zero_hit;

err4 = sem_minus2_miss;
err5 = sem_minus1_miss;
err6 = sem_zero_miss;


plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
hold on
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])

xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'hit','','miss'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Delta IBI')
title('ACD (hit vs miss): Day 1')
set(gca,'FontSize',36)
hold off
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'acd_min2_hit','acd_min1_hit','acd_zero_hit','acd_min2_miss','acd_min1_miss','acd_zero_miss'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Detection','Time'});
withinDesign.Detection = categorical(withinDesign.Detection);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'acd_min2_hit-acd_zero_miss ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Detection*Time');
disp(anovaTable(AT, 'Value'));

%% ACD progression plots -- High vs. Low confidence
nSubj=19;
for s=1:nSubj
    conf = setdiff(1:23,[4,13,21,23]);
    mean_acd_s2_highC(s,1) = mean(acd(conf(s)).acd_s2_highC,'omitnan'); mean_acd_s2_lowC(s,1) = mean(acd(conf(s)).acd_s2_lowC,'omitnan');
    mean_acd_s1_highC(s,1) = mean(acd(conf(s)).acd_s1_highC,'omitnan'); mean_acd_s1_lowC(s,1) = mean(acd(conf(s)).acd_s1_lowC,'omitnan');
    mean_acd_s0_highC(s,1) = mean(acd(conf(s)).acd_s0_highC,'omitnan'); mean_acd_s0_lowC(s,1) = mean(acd(conf(s)).acd_s0_lowC,'omitnan');
%     mean_acd_sp1_highC(s,1) = mean(acd(s).acd_s_p1_highC,'omitnan'); mean_acd_sp1_lowC(s,1) = mean(acd(s).acd_s_p1_lowC,'omitnan');
    
    sem_acd_s2_highC(s,1) = std(acd(conf(s)).acd_s2_highC,'omitnan')/sqrt(size(acd(conf(s)).acd_s2_highC,2)); sem_acd_s2_lowC(s,1) = std(acd(conf(s)).acd_s2_lowC,'omitnan')/sqrt(size(acd(conf(s)).acd_s2_lowC,2));
    sem_acd_s1_highC(s,1) = std(acd(conf(s)).acd_s1_highC,'omitnan')/sqrt(size(acd(conf(s)).acd_s1_highC,2)); sem_acd_s1_lowC(s,1) = std(acd(conf(s)).acd_s1_lowC,'omitnan')/sqrt(size(acd(conf(s)).acd_s1_lowC,2));
    sem_acd_s0_highC(s,1) = std(acd(conf(s)).acd_s0_highC,'omitnan')/sqrt(size(acd(conf(s)).acd_s0_highC,2)); sem_acd_s0_lowC(s,1) = std(acd(conf(s)).acd_s0_lowC,'omitnan')/sqrt(size(acd(conf(s)).acd_s0_lowC,2));
%     sem_acd_sp1_highC(s,1) = std(acd(s).acd_s_p1_highC,'omitnan')/sqrt(size(acd(s).acd_s_p1_highC,2)); sem_acd_sp1_lowC(s,1) = std(acd(s).acd_s_p1_lowC,'omitnan')/sqrt(size(acd(s).acd_s_p1_lowC,2));
    
end
var11 = mean_acd_s2_highC;
var12 = mean_acd_s1_highC;
var13 = mean_acd_s0_highC;
var21 = mean_acd_s2_lowC;
var22 = mean_acd_s1_lowC;
var23 = mean_acd_s0_lowC;
var1 = mean([var11 var12 var13],2);
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];
var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_highC = std(var11)/sqrt(length(var11)); sem_minus2_lowC = std(var21)/sqrt(length(var21));
sem_minus1_highC = std(var12)/sqrt(length(var12)); sem_minus1_lowC = std(var22)/sqrt(length(var22));
sem_zero_highC = std(var13)/sqrt(length(var13)); sem_zero_lowC = std(var23)/sqrt(length(var23));

%%% plot the IBI progression over time in Fast and Slow responses separately
figure
times = [1:3];

err1 = sem_minus2_highC;
err2 = sem_minus1_highC;
err3 = sem_zero_highC;

err4 = sem_minus2_lowC;
err5 = sem_minus1_lowC;
err6 = sem_zero_lowC;


plot(times, [mean(var11,'omitnan') mean(var12,'omitnan') mean(var13,'omitnan')],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11,'omitnan') mean(var12,'omitnan') mean(var13,'omitnan')],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
hold on
plot(times, [mean(var21,'omitnan') mean(var22,'omitnan') mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])

xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'ACD-highC','','ACD-lowC'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Delta IBI')
title('ACD (highC vs lowC): Day 1')
set(gca,'FontSize',36)
hold off
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'acd_min2_highC','acd_min1_highC','acd_zero_highC','acd_min2_lowC','acd_min1_lowC','acd_zero_lowC'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Confidence','Time'});
withinDesign.Confidence = categorical(withinDesign.Confidence);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'acd_min2_highC-acd_zero_lowC ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Confidence*Time');
disp(anovaTable(AT, 'Value'));

%% s2 --> s1 --> s0 ACD one way progress
var11 = mean_acd_s2; var12 = mean_acd_s1; var13 = mean_acd_s0;
var1 = mean([var11 var12 var13],2);
var11 = var11(~isoutlier(var1)); var12 = var12(~isoutlier(var1)); var13 = var13(~isoutlier(var1));

%%% calculate group level SEMs
sem_minus2_acd = std(var11)/sqrt(length(var11));
sem_minus1_acd = std(var12)/sqrt(length(var12));
sem_zero_acd = std(var13)/sqrt(length(var13));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_acd);
err2 = mean(sem_minus1_acd);
err3 = mean(sem_zero_acd);

% Time points
times = [1 2 3];
%%%
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0 0])

xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);

xticklabels({'Stim-2','Stim-1','Stim'})
legend({'ACD'},'Location','northwest')
xlabel('time')
ylabel('ACD (ms)')
title('ACD in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

% Repeated measures model
T = table(var11, var12, var13, 'VariableNames', {'s2', 's1', 's0'});
rm = fitrm(T, 's2-s0 ~ 1', 'WithinDesign', table([1 2 3]','VariableNames',{'Condition'}));

% Run repeated measures ANOVA
ranovatbl = ranova(rm)
mauchly(rm)
% View results
disp(ranovatbl)

%% s2 --> s1 --> s0 progression of phi, A, B

%% PHI
var11 = [rsa_ibi_modelData_BBB.phi_s2]'; var12 = [rsa_ibi_modelData_BBB.phi_s1]'; var13 = [rsa_ibi_modelData_BBB.phi_s0]';
var1 = mean([var11 var12 var13],2);
var11 = var11(~isoutlier(var1)); var12 = var12(~isoutlier(var1)); var13 = var13(~isoutlier(var1));
allVar = [var11;var12;var13];
% var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
% var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
% var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));

var10 = [baseRSA.phi]';
%%% calculate group level SEMs
sem_base_phi = std(var10)/sqrt(length(var10));
sem_minus2_phi = std(var11)/sqrt(length(var11));
sem_minus1_phi = std(var12)/sqrt(length(var12));
sem_zero_phi = std(var13)/sqrt(length(var13));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:4];

err1 = mean(sem_base_phi);
err2 = mean(sem_minus2_phi);
err3 = mean(sem_minus1_phi);
err4 = mean(sem_zero_phi);

% Time points
times = [1 2 3 4];
%%%
plot(times, [mean(var10) mean(var11) mean(var12) mean(var13)],'Color',[0 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var10) mean(var11) mean(var12) mean(var13)],[err1,err2,err3,err4],'-o', 'LineWidth', 2, 'Color', [0 0 0])

xlim([0.5 4.5])
% ylim([840 950])
xticks([1 2 3 4])
% yticks([860 880 900 920 940]);

xticklabels({'Rest','Stim-2','Stim-1','Stim'})
legend({'phi'},'Location','northwest')
xlabel('time')
ylabel('Sine Fits: Phase Shift')
title('Phase Shift in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

% Repeated measures model
T = table(var10, var11, var12, var13, 'VariableNames', {'rest','s2', 's1', 's0'});
rm = fitrm(T, 'rest-s0 ~ 1', 'WithinDesign', table([1 2 3 4]','VariableNames',{'Condition'}));

% Run repeated measures ANOVA
ranovatbl = ranova(rm)

% View results
disp(ranovatbl)

%% Amplitude
var11 = [rsa_ibi_modelData_BBB.A_s2]'; var12 = [rsa_ibi_modelData_BBB.A_s1]'; var13 = [rsa_ibi_modelData_BBB.A_s0]';
var1 = mean([var11 var12 var13],2);
var11 = var11(~isoutlier(var1)); var12 = var12(~isoutlier(var1)); var13 = var13(~isoutlier(var1));
% allVar = [var11;var12;var13];
% var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
% var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
% var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));

var10 = [baseRSA.A]';
%%% calculate group level SEMs
sem_base_A = std(var10)/sqrt(length(var10));
sem_minus2_A = std(var11)/sqrt(length(var11));
sem_minus1_A = std(var12)/sqrt(length(var12));
sem_zero_A = std(var13)/sqrt(length(var13));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
subplot(1,2,1)
times = [1:4];

err1 = mean(sem_base_A);
err2 = mean(sem_minus2_A);
err3 = mean(sem_minus1_A);
err4 = mean(sem_zero_A);

% Time points
times = [1 2 3 4];
%%%
plot(times, [mean(var10) mean(var11) mean(var12) mean(var13)],'Color',[0 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var10) mean(var11) mean(var12) mean(var13)],[err1,err2,err3,err4],'-o', 'LineWidth', 2, 'Color', [0 0 0])

xlim([0.5 4.5])
% ylim([840 950])
xticks([1 2 3 4])
% yticks([860 880 900 920 940]);

xticklabels({'Rest','Stim-2','Stim-1','Stim'})
legend({'A'},'Location','northwest')
% xlabel('time')
ylabel('Sine Fits: Amplitude')
title('Amplitude in Stimulus Anticipation')
set(gca,'FontSize',36)
% hold off

% Repeated measures model
T = table(var10, var11, var12, var13, 'VariableNames', {'rest','s2', 's1', 's0'});
rm = fitrm(T, 'rest-s0 ~ 1', 'WithinDesign', table([1 2 3 4]','VariableNames',{'Condition'}));

% Run repeated measures ANOVA
ranovatbl = ranova(rm)

% View results
disp(ranovatbl)

%%%%%%%%%%%%%
%%% B (vertical offset)
%%%%%%%%%%%%%

var11 = [rsa_ibi_modelData_BBB.B_s2]'; var12 = [rsa_ibi_modelData_BBB.B_s1]'; var13 = [rsa_ibi_modelData_BBB.B_s0]';
var1 = mean([var11 var12 var13],2);
var11 = var11(~isoutlier(var1)); var12 = var12(~isoutlier(var1)); var13 = var13(~isoutlier(var1));
% allVar = [var11;var12;var13];
% var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
% var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
% var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));

var10 = [baseRSA.B]';
%%% calculate group level SEMs
sem_base_B = std(var10)/sqrt(length(var10));
sem_minus2_B = std(var11)/sqrt(length(var11));
sem_minus1_B = std(var12)/sqrt(length(var12));
sem_zero_B = std(var13)/sqrt(length(var13));

%%% plot the vertical offset progression over time in Hit vs Miss trials
subplot(1,2,2)
times = [1:4];

err1 = mean(sem_base_B);
err2 = mean(sem_minus2_B);
err3 = mean(sem_minus1_B);
err4 = mean(sem_zero_B);

% Time points
times = [1 2 3 4];
%%%
plot(times, [mean(var10) mean(var11) mean(var12) mean(var13)],'Color',[0 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var10) mean(var11) mean(var12) mean(var13)],[err1,err2,err3,err4],'-o', 'LineWidth', 2, 'Color', [0 0 0])

xlim([0.5 4.5])
% ylim([840 950])
xticks([1 2 3 4])
% yticks([860 880 900 920 940]);

xticklabels({'Rest','Stim-2','Stim-1','Stim'})
legend({'B'},'Location','northwest')
% xlabel('time')
ylabel('Sine Fits: Vertical Offset')
title('Vertical Offset in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

% Repeated measures model
T = table(var10, var11, var12, var13, 'VariableNames', {'rest','s2', 's1', 's0'});
rm = fitrm(T, 'rest-s0 ~ 1', 'WithinDesign', table([1 2 3 4]','VariableNames',{'Condition'}));

% Run repeated measures ANOVA
ranovatbl = ranova(rm)

% View results
disp(ranovatbl)

%% Sine fit (phi -- phase-shift) DETECTION
var11 = [rsa_ibi_modelData_BBB.phi_hit_s2]';
var12 = [rsa_ibi_modelData_BBB.phi_hit_s1]';
var13 = [rsa_ibi_modelData_BBB.phi_hit_s0]';
var1 = mean([var11 var12 var13],2);

var21 = [rsa_ibi_modelData_BBB.phi_miss_s2]';
var22 = [rsa_ibi_modelData_BBB.phi_miss_s1]';
var23 = [rsa_ibi_modelData_BBB.phi_miss_s0]';
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];

var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_phiHit = std(var11)/sqrt(length(var11)); sem_minus2_phiMiss = std(var21)/sqrt(length(var21));
sem_minus1_phiHit = std(var12)/sqrt(length(var12)); sem_minus1_phiMiss = std(var22)/sqrt(length(var22));
sem_zero_phiHit = std(var13)/sqrt(length(var13)); sem_zero_phiMiss = std(var23)/sqrt(length(var23));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_phiHit);
err2 = mean(sem_minus1_phiHit);
err3 = mean(sem_zero_phiHit);

err4 = mean(sem_minus2_phiMiss);
err5 = mean(sem_minus1_phiMiss);
err6 = mean(sem_zero_phiMiss);

% Time points
times = [1 2 3];
%%% Hits
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
%%% Misses
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])
xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'phi-hit','','phi-miss'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Sine Fits: Phase Shift')
title('Phase Shift in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'phi_min2_hit','phi_min1_hit','phi_zero_hit','phi_min2_miss','phi_min1_miss','phi_zero_miss'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Detection','Time'});
withinDesign.Detection = categorical(withinDesign.Detection);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'phi_min2_hit-phi_zero_miss ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Detection*Time');
disp(anovaTable(AT, 'Value'));

%% Sine fit (phi -- phase-shift) CONFIDENCE
var11 = [rsa_ibi_modelData_BBB.phi_highC_s2]'; var11 = var11(~isnan(var11));
var12 = [rsa_ibi_modelData_BBB.phi_highC_s1]'; var12 = var12(~isnan(var12));
var13 = [rsa_ibi_modelData_BBB.phi_highC_s0]'; var13 = var13(~isnan(var13));
var1 = mean([var11 var12 var13],2);

var21 = [rsa_ibi_modelData_BBB.phi_lowC_s2]'; var21 = var21(~isnan(var21));
var22 = [rsa_ibi_modelData_BBB.phi_lowC_s1]'; var22 = var22(~isnan(var22));
var23 = [rsa_ibi_modelData_BBB.phi_lowC_s0]'; var23 = var23(~isnan(var23));
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];

var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_phiHit = std(var11)/sqrt(length(var11)); sem_minus2_phiMiss = std(var21)/sqrt(length(var21));
sem_minus1_phiHit = std(var12)/sqrt(length(var12)); sem_minus1_phiMiss = std(var22)/sqrt(length(var22));
sem_zero_phiHit = std(var13)/sqrt(length(var13)); sem_zero_phiMiss = std(var23)/sqrt(length(var23));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_phiHit);
err2 = mean(sem_minus1_phiHit);
err3 = mean(sem_zero_phiHit);

err4 = mean(sem_minus2_phiMiss);
err5 = mean(sem_minus1_phiMiss);
err6 = mean(sem_zero_phiMiss);

% Time points
times = [1 2 3];
%%% Hits
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
%%% Misses
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])
xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'phi-highC','','phi-lowC'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Sine Fits: Phase Shift')
title('Phase Shift in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'phi_min2_highC','phi_min1_highC','phi_zero_highC','phi_min2_lowC','phi_min1_lowC','phi_zero_lowC'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Confidence','Time'});
withinDesign.Confidence = categorical(withinDesign.Confidence);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'phi_min2_highC-phi_zero_lowC ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Confidence*Time');
disp(anovaTable(AT, 'Value'));

%% Sine fit (A -- amplitude) progression in Hit vs Miss trials over the course of stimulus anticipation
var11 = [rsa_ibi_modelData_BBB.A_hit_s2]';
var12 = [rsa_ibi_modelData_BBB.A_hit_s1]';
var13 = [rsa_ibi_modelData_BBB.A_hit_s0]';
var1 = mean([var11 var12 var13],2);

var21 = [rsa_ibi_modelData_BBB.A_miss_s2]';
var22 = [rsa_ibi_modelData_BBB.A_miss_s1]';
var23 = [rsa_ibi_modelData_BBB.A_miss_s0]';
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];

var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_AHit = std(var11)/sqrt(length(var11)); sem_minus2_AMiss = std(var21)/sqrt(length(var21));
sem_minus1_AHit = std(var12)/sqrt(length(var12)); sem_minus1_AMiss = std(var22)/sqrt(length(var22));
sem_zero_AHit = std(var13)/sqrt(length(var13)); sem_zero_AMiss = std(var23)/sqrt(length(var23));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_AHit);
err2 = mean(sem_minus1_AHit);
err3 = mean(sem_zero_AHit);

err4 = mean(sem_minus2_AMiss);
err5 = mean(sem_minus1_AMiss);
err6 = mean(sem_zero_AMiss);

% Time points
times = [1 2 3];
%%% Hits
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
%%% Misses
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])
xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'A-hit','','A-miss'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Sine Fits: Amplitude Shift')
title('Amplitude Shift in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'A_min2_hit','A_min1_hit','A_zero_hit','A_min2_miss','A_min1_miss','A_zero_miss'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Detection','Time'});
withinDesign.Detection = categorical(withinDesign.Detection);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'A_min2_hit-A_zero_miss ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Detection*Time');
disp(anovaTable(AT, 'Value'));

%%% Friedman
% Detection effect at time1
p_det_1 = signrank(var11, var21)  % Wilcoxon signed-rank test

% Detection effect at time2
p_det_2 = signrank(var12, var22)

% Detection effect at time3
p_det_3 = signrank(var13, var23)

%% Pre-Stimulus Amplitude Comparison: Hit vs Miss: signrank([rsa_ibi_modelData.A_hit]',[rsa_ibi_modelData.A_miss]')
%%% 
hit_A = [rsa_ibi_modelData.A_hit]';
miss_A = [rsa_ibi_modelData.A_miss]';
var1=hit_A;
var2=miss_A;

var12 = [{var1}; {var2}];
[p h] = signrank(var1,var2); %ttest also gives a very similar result

figure
h1 = rm_raincloud(var12,[0 0.8 0],0,'ks',[],0.01);

h1.p{2,1}.FaceColor=[0.8 0 0]; %face color for the second patch
h1.s{2,1}.MarkerFaceColor=[0.8 0 0]; %marker taskColor for the second category
h1.l(1,1).Visible="off";
h1.m(2,1).MarkerFaceColor = [0.8 0 0];

yticklabels({'Miss','Hit'}) %inverted because of the rm_raincloud function plot orientation!

% Also match individual data points
hold on
for i=1:size(var1,1)
    X = [h1.s{1,1}.XData(1,i) h1.s{2,1}.XData(1,i)];
    Y = [h1.s{1,1}.YData(1,i) h1.s{2,1}.YData(1,i)];
    plot(X,Y,"k-");
end

% Plot three lines as your significance line
% 1
hold on
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5];
Yv = [mean(h1.s{1,1}.YData) mean(h1.s{2,1}.YData)];
plot(Xv,Yv,"k-");
hold on
% 2
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+3 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5];
Yv = [mean(h1.s{1,1}.YData) mean(h1.s{1,1}.YData)];
plot(Xv,Yv,"k-");
% 3
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+3 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5];
Yv = [mean(h1.s{2,1}.YData) mean(h1.s{2,1}.YData)];
plot(Xv,Yv,"k-");

hold on
pF = round(p,2);
text(max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+10,(mean(h1.s{1,1}.YData)+mean(h1.s{2,1}.YData))/2,['p = ',num2str(pF)],"FontSize",36,"HorizontalAlignment","center",'FontWeight','bold');
% xlim([0 80])
% xticks([-40 -30 -20 -10 0 10 20 30 40 50 60])
set(gca,'FontSize',36,'FontWeight','Bold')
set(gcf, 'Position', get(0, 'Screensize'));
xlabel('RSA (A)') %we need to invert the x and y for all the labels and ticks as well due to rainplots being scripted in a rotated way!
title(['RSA Amplitude in Stimulus Anticipation: Hit vs Miss'])
hold off
%% Sine fit (A -- amplitude) CONFIDENCE
var11 = [rsa_ibi_modelData_BBB.A_highC_s2]'; var11 = var11(~isnan(var11));
var12 = [rsa_ibi_modelData_BBB.A_highC_s1]'; var12 = var12(~isnan(var12));
var13 = [rsa_ibi_modelData_BBB.A_highC_s0]'; var13 = var13(~isnan(var13));
var1 = mean([var11 var12 var13],2);

var21 = [rsa_ibi_modelData_BBB.A_lowC_s2]'; var21 = var21(~isnan(var21));
var22 = [rsa_ibi_modelData_BBB.A_lowC_s1]'; var22 = var22(~isnan(var22));
var23 = [rsa_ibi_modelData_BBB.A_lowC_s0]'; var23 = var23(~isnan(var23));
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];

var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_AHit = std(var11)/sqrt(length(var11)); sem_minus2_AMiss = std(var21)/sqrt(length(var21));
sem_minus1_AHit = std(var12)/sqrt(length(var12)); sem_minus1_AMiss = std(var22)/sqrt(length(var22));
sem_zero_AHit = std(var13)/sqrt(length(var13)); sem_zero_AMiss = std(var23)/sqrt(length(var23));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_AHit);
err2 = mean(sem_minus1_AHit);
err3 = mean(sem_zero_AHit);

err4 = mean(sem_minus2_AMiss);
err5 = mean(sem_minus1_AMiss);
err6 = mean(sem_zero_AMiss);

% Time points
times = [1 2 3];
%%% Hits
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
%%% Misses
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])
xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'A-highC','','A-lowC'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Sine Fits: Amplitude')
title('Amplitude in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'A_min2_highC','A_min1_highC','A_zero_highC','A_min2_lowC','A_min1_lowC','A_zero_lowC'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Confidence','Time'});
withinDesign.Confidence = categorical(withinDesign.Confidence);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'A_min2_highC-A_zero_lowC ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Confidence*Time');
disp(anovaTable(AT, 'Value'));

%%% Friedman
highC = [var11;var12;var13];
lowC = [var21;var22;var23];
tableDetect = [highC lowC];
[p, tbl, stats] = friedman(tableDetect,size(var11,1))

%% Sine fit (B -- offset) progression in Hit vs Miss trials over the course of stimulus anticipation
var11 = [rsa_ibi_modelData_BBB.B_hit_s2]'-[rsa_ibi_modelData_BBB.B_hit_s2]';
var12 = [rsa_ibi_modelData_BBB.B_hit_s1]'-[rsa_ibi_modelData_BBB.B_hit_s2]';
var13 = [rsa_ibi_modelData_BBB.B_hit_s0]'-[rsa_ibi_modelData_BBB.B_hit_s2]';
var1 = mean([var11 var12 var13],2);

var21 = [rsa_ibi_modelData_BBB.B_miss_s2]'-[rsa_ibi_modelData_BBB.B_miss_s2]';
var22 = [rsa_ibi_modelData_BBB.B_miss_s1]'-[rsa_ibi_modelData_BBB.B_miss_s2]';
var23 = [rsa_ibi_modelData_BBB.B_miss_s0]'-[rsa_ibi_modelData_BBB.B_miss_s2]';
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];

% var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
% var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

% var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
% var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
% var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
% var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
% var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
% var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_BHit = std(var11)/sqrt(length(var11)); sem_minus2_BMiss = std(var21)/sqrt(length(var21));
sem_minus1_BHit = std(var12)/sqrt(length(var12)); sem_minus1_BMiss = std(var22)/sqrt(length(var22));
sem_zero_BHit = std(var13)/sqrt(length(var13)); sem_zero_BMiss = std(var23)/sqrt(length(var23));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_BHit);
err2 = mean(sem_minus1_BHit);
err3 = mean(sem_zero_BHit);

err4 = mean(sem_minus2_BMiss);
err5 = mean(sem_minus1_BMiss);
err6 = mean(sem_zero_BMiss);

% Time points
times = [1 2 3];
%%% Hits
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
%%% Misses
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])
xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'B-hit','','B-miss'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Sine Fits: Vertical Offset')
title('Vertical Offset in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'B_min2_hit','B_min1_hit','B_zero_hit','B_min2_miss','B_min1_miss','B_zero_miss'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Detection','Time'});
withinDesign.Detection = categorical(withinDesign.Detection);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'B_min2_hit-B_zero_miss ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Detection*Time');
disp(anovaTable(AT, 'Value'));

%% Sine fit (B -- offset) in Hit vs Miss trials over the course of stimulus anticipation: s2 vs s0 (total drift)
hit_dOffset = var13-var11;
miss_dOffset = var23-var21;
var1=hit_dOffset;
var2=miss_dOffset;

var12 = [{var1}; {var2}];
[p h stats] = signrank(var1,var2); %ttest also gives a very similar result

figure
h1 = rm_raincloud(var12,[0 0.8 0],0,'ks',[],0.01);

h1.p{2,1}.FaceColor=[0.8 0 0]; %face color for the second patch
h1.s{2,1}.MarkerFaceColor=[0.8 0 0]; %marker taskColor for the second category
h1.l(1,1).Visible="off";
h1.m(2,1).MarkerFaceColor = [0.8 0 0];

yticklabels({'Miss','Hit'}) %inverted because of the rm_raincloud function plot orientation!

% Also match individual data points
hold on
for i=1:size(var1,1)
    X = [h1.s{1,1}.XData(1,i) h1.s{2,1}.XData(1,i)];
    Y = [h1.s{1,1}.YData(1,i) h1.s{2,1}.YData(1,i)];
    plot(X,Y,"k-");
end

% Plot three lines as your significance line
% 1
hold on
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1];
Yv = [mean(h1.s{1,1}.YData) mean(h1.s{2,1}.YData)];
plot(Xv,Yv,"k-");
hold on
% 2
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.06 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1];
Yv = [mean(h1.s{1,1}.YData) mean(h1.s{1,1}.YData)];
plot(Xv,Yv,"k-");
% 3
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.06 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1];
Yv = [mean(h1.s{2,1}.YData) mean(h1.s{2,1}.YData)];
plot(Xv,Yv,"k-");

hold on
pF = round(p,2);
text(max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2,(mean(h1.s{1,1}.YData)+mean(h1.s{2,1}.YData))/2,['p = ',num2str(pF)],"FontSize",36,"HorizontalAlignment","center",'FontWeight','bold');
% xlim([-40 60])
% xticks([-40 -30 -20 -10 0 10 20 30 40 50 60])
set(gca,'FontSize',36,'FontWeight','Bold')
set(gcf, 'Position', get(0, 'Screensize'));
xlabel('Delta Offset') %we need to invert the x and y for all the labels and ticks as well due to rainplots being scripted in a rotated way!
title(['Non-Respiratory Cardiac Deceleration: Hit vs Miss'])
hold off
%%
figure
subplot(1,2,1)
h1=daboxplot([{[rsa_ibi_modelData.B_hit]'} {[rsa_ibi_modelData.B_miss]'}],'scatter',1,'scattersize',300,'mean',1,'colors',[0 0.6 0; 0.6 0 0]);
% Also match individual data points
hold on
for i=1:size(var1,1)
    X = [h1.sc(1,1).XData(1,i) h1.sc(1,2).XData(1,i)];
    Y = [h1.sc(1,1).YData(1,i) h1.sc(1,2).YData(1,i)];
    plot(X,Y,"k--",'LineWidth',0.5);
end
% ylim([-32 50])
% yticks([-0.2 0 0.2 0.4 0.6])
xticklabels({'Hit','Miss'})
ylabel('delta B (norm)')
set(gca,'FontSize',36)
% set(gcf, 'Position', get(0, 'Screensize'));
[p h stats] = signrank([rsa_ibi_modelData.B_hit]',[rsa_ibi_modelData.B_miss]')

subplot(1,2,2)
h1=daboxplot([{hit_dOffset} {miss_dOffset}],'scatter',1,'scattersize',300,'mean',1,'colors',[0 0.6 0; 0.6 0 0]);
% Also match individual data points
hold on
for i=1:size(var1,1)
    X = [h1.sc(1,1).XData(1,i) h1.sc(1,2).XData(1,i)];
    Y = [h1.sc(1,1).YData(1,i) h1.sc(1,2).YData(1,i)];
    plot(X,Y,"k--",'LineWidth',0.5);
end
% ylim([-32 50])
% yticks([-0.2 0 0.2 0.4 0.6])
xticklabels({'Hit','Miss'})
ylabel('delta B (norm)')
set(gca,'FontSize',36)
% set(gcf, 'Position', get(0, 'Screensize'));
% [p h stats] = signrank([rsa_ibi_modelData.B_hit]',[rsa_ibi_modelData.B_miss]')

%% Diff between hit and miss
var1=hit_dOffset;
% var1(11) = [];
var2=miss_dOffset;
% var2(11) = [];
subplot(1,2,2)
raincloud_plot(var1 - var2,'box_on', 1,'color',[0.5 0.5 0.5]);
hold on
l1 = xline(0,'color',[0.2 0.2 0.2],'LineWidth',3,'LineStyle','--','Label','Origin','LabelHorizontalAlignment','left','LabelOrientation','horizontal','FontSize',32);
% l2 = xline(median(hit_avgPrePupil_mm - miss_avgPrePupil_mm),'color',[0.8 0 0],'LineWidth',3,'LineStyle','--');
yticks([]);
xlabel('Delta Offset Difference');
% xlim([-20 40])
title(['Hit Minus Miss: Vertical Offset Increase'])
set(gca,'FontSize',36)
% xticks([-20 -10 0 10 20 30 40])

%% Sine fit (B -- vertical offset) CONFIDENCE
var11 = [rsa_ibi_modelData_BBB.B_highC_s2]'-[rsa_ibi_modelData_BBB.B_highC_s2]'; var11 = var11(~isnan(var11));
var12 = [rsa_ibi_modelData_BBB.B_highC_s1]'-[rsa_ibi_modelData_BBB.B_highC_s2]'; var12 = var12(~isnan(var12));
var13 = [rsa_ibi_modelData_BBB.B_highC_s0]'-[rsa_ibi_modelData_BBB.B_highC_s2]'; var13 = var13(~isnan(var13));
var1 = mean([var11 var12 var13],2);

var21 = [rsa_ibi_modelData_BBB.B_lowC_s2]'-[rsa_ibi_modelData_BBB.B_lowC_s2]'; var21 = var21(~isnan(var21));
var22 = [rsa_ibi_modelData_BBB.B_lowC_s1]'-[rsa_ibi_modelData_BBB.B_lowC_s2]'; var22 = var22(~isnan(var22));
var23 = [rsa_ibi_modelData_BBB.B_lowC_s0]'-[rsa_ibi_modelData_BBB.B_lowC_s2]'; var23 = var23(~isnan(var23));
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];
% 
var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

% var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
% var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
% var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
% var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
% var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
% var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

%%% calculate group level SEMs
sem_minus2_BHit = std(var11)/sqrt(length(var11)); sem_minus2_BMiss = std(var21)/sqrt(length(var21));
sem_minus1_BHit = std(var12)/sqrt(length(var12)); sem_minus1_BMiss = std(var22)/sqrt(length(var22));
sem_zero_BHit = std(var13)/sqrt(length(var13)); sem_zero_BMiss = std(var23)/sqrt(length(var23));

%%% plot the vertical offset progression over time in Hit vs Miss trials
figure
times = [1:3];

err1 = mean(sem_minus2_BHit);
err2 = mean(sem_minus1_BHit);
err3 = mean(sem_zero_BHit);

err4 = mean(sem_minus2_BMiss);
err5 = mean(sem_minus1_BMiss);
err6 = mean(sem_zero_BMiss);

% Time points
times = [1 2 3];
%%% Hits
plot(times, [mean(var11) mean(var12) mean(var13)],'Color',[0 0.6 0],'LineWidth',4);
hold on
errorbar(times,[mean(var11) mean(var12) mean(var13)],[err1,err2,err3],'-o', 'LineWidth', 2, 'Color', [0 0.6 0])
%%% Misses
plot(times, [mean(var21) mean(var22) mean(var23)],'Color',[0.6 0 0],'LineWidth',4);
hold on
errorbar(times,[mean(var21) mean(var22) mean(var23)],[err4,err5,err6],'-o', 'LineWidth', 2, 'Color', [0.6 0 0])
xlim([0.5 3.5])
% ylim([840 950])
xticks([1 2 3])
% yticks([860 880 900 920 940]);


xticklabels({'Stim-2','Stim-1','Stim'})
legend({'B-highC','','B-lowC'},'Location','northwest')
% xlabel('RR-interval timings')
ylabel('Sine Fits: Vertical Offset')
title('Vertical Offset in Stimulus Anticipation')
set(gca,'FontSize',36)
hold off

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts

%%% two-way repeated measures ANOVA on the interaction between vertical offsets on the sine model AND stimulus anticipatory period!
% 
intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'B_min2_highC','B_min1_highC','B_zero_highC','B_min2_lowC','B_min1_lowC','B_zero_lowC'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'Confidence','Time'});
withinDesign.Confidence = categorical(withinDesign.Confidence);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'B_min2_highC-B_zero_lowC ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'Confidence*Time');
disp(anovaTable(AT, 'Value'));

%% exhale vs inhale (relevant "breathingData" dataset is "acd_exhalationVsInhalation.mat")

acd_exh_s2=[breathingData.mean_acd_exhale_s2]'; acd_exh_s1= [breathingData.mean_acd_exhale_s1]'; acd_exh_s0 = [breathingData.mean_acd_exhale_s0]';
preStim_acd_exh_mean = mean([acd_exh_s2 acd_exh_s1 acd_exh_s0],2);
acd_inh_s2=[breathingData.mean_acd_inhale_s2]'; acd_inh_s1= [breathingData.mean_acd_inhale_s1]'; acd_inh_s0 = [breathingData.mean_acd_inhale_s0]';
preStim_acd_inh_mean = mean([acd_inh_s2 acd_inh_s1 acd_inh_s0],2);

preStim_RSA = preStim_acd_exh_mean - preStim_acd_inh_mean;
preStim_nonRSA = preStim_acd_exh_mean + preStim_acd_inh_mean;

breathingMatrix = [[breathingData.mean_acd_exhale_s2]' [breathingData.mean_acd_inhale_s2]' [breathingData.mean_acd_exhale_s1]' [breathingData.mean_acd_inhale_s1]'...
    [breathingData.mean_acd_exhale_s0]' [breathingData.mean_acd_inhale_s0]'];
breathingArray = mat2cell(breathingMatrix, size(breathingMatrix, 1), ones(1, size(breathingMatrix, 2)));
groups = [1 2 1 2 1 2];
colors = [0 0.5 0.5; 0.5 0.5 0.5; 0 0.5 0.5; 0.5 0.5 0.5; 0 0.5 0.5; 0.5 0.5 0.5];
xtlab = {'exh-2','inh-2','exh-1','inh-1','exh-0','inh-0'};

figure
dabarplot(breathingArray,'groups',groups,'colors',colors,'scatter',1,'scattersize',40,'xtlabels',xtlab);
ylabel('Anticipatory Cardiac Deceleration')
xlabel('Exhalation (exh) vs. Inhalation (inh)')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
set(gca, 'FontSize', 36);

var11 = [breathingData.mean_acd_exhale_s2]'; var12 = [breathingData.mean_acd_exhale_s1]'; var13 = [breathingData.mean_acd_exhale_s0]';
var1 = mean([var11 var12 var13],2);
var21 = [breathingData.mean_acd_inhale_s2]'; var22 = [breathingData.mean_acd_inhale_s1]'; var23 = [breathingData.mean_acd_inhale_s0]';
var2 = mean([var21 var22 var23],2);
allVar = [var11;var12;var13;var21;var22;var23];
var11 = var11(~isoutlier(var1) & ~isoutlier(var2)); var12 = var12(~isoutlier(var1) & ~isoutlier(var2)); var13 = var13(~isoutlier(var1) & ~isoutlier(var2));
var21 = var21(~isoutlier(var1) & ~isoutlier(var2)); var22 = var22(~isoutlier(var1) & ~isoutlier(var2)); var23 = var23(~isoutlier(var1) & ~isoutlier(var2));

%%%normalize
var13 = (var13 - min(allVar)) / (max(allVar)-min(allVar));
var23 = (var23 - min(allVar)) / (max(allVar)-min(allVar));
var12 = (var12 - min(allVar)) / (max(allVar)-min(allVar));
var22 = (var22 - min(allVar)) / (max(allVar)-min(allVar));
var11 = (var11 - min(allVar)) / (max(allVar)-min(allVar));
var21 = (var21 - min(allVar)) / (max(allVar)-min(allVar));

intTable = table(var11, var12, var13, var21, var22, var23);
% 
intTable.Properties.VariableNames = {'exh_s2','exh_s1','exh_s0','inh_s2','inh_s1','inh_s0'};
% 
% % create the within-subjects design
withinDesign = table([1 1 1 2 2 2]', [1 2 3 1 2 3]', 'VariableNames',{'BrPhase','Time'});
withinDesign.BrPhase = categorical(withinDesign.BrPhase);
withinDesign.Time = categorical(withinDesign.Time);
% 
% %create repeated measures model
rm = fitrm(intTable, 'exh_s2-inh_s0 ~ 1', 'WithinDesign', withinDesign);
% 
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK
AT = ranova(rm, 'WithinModel', 'BrPhase*Time');
disp(anovaTable(AT, 'Value'));

%% multiple linear regression computations
X = [ones(size([rsa_ibi_modelData.A],2),1) [rsa_ibi_modelData.A]' [baseRSA.A]' [rsa_ibi_modelData.A]'.*[baseRSA.A]']; %baseline and anticipatory RSA
y = (mean_acd_s0+mean_acd_s1+mean_acd_s2)/3; %mean ACD
[b,bint,r,rint,stats] = regress(y,X)

%% Linear regression computations (you can modify the x and y)
x = [baseRSA.A]';
% x = x(setdiff(1:23,find(isoutlier(baselineHFs_brSubj))));
y = [rsa_ibi_modelData.B]';  %(setdiff(1:23,find(isoutlier(baselineHFs_brSubj))));

x_c = x(~isoutlier(x) & ~isoutlier(y));
y_c = y(~isoutlier(x) & ~isoutlier(y));

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
plot(x_c, y_fit, '-','LineWidth',10);     % Plot the fitted line
hold off;
xlabel('Resting-State RSA');
ylabel('Anticipatory B');
title('RSA vs. B');
legend('Data', 'Fitted Line');

% Add text to the northwest corner
text_loc_x = min(x_c); % Use the minimum x value for positioning
text_loc_y = max(y_c); % Use the maximum y value for positioning

[r, p_value] = corr(x_c, y_c)

% Calculate correlation coefficient and p-value
% Fit robust linear models
mdl_huber = fitlm(x_c, y_c, 'RobustOpts', 'huber')
cor = sqrt(mdl_huber.Rsquared.Ordinary)

set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
text(text_loc_x, text_loc_y, ...
    sprintf('slope = %.2f\np = %.4f', mdl_huber.Coefficients.Estimate(2), mdl_huber.Coefficients.pValue(2)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline', ...
    'FontSize', 48, 'BackgroundColor', 'white', 'EdgeColor', 'black');
set(gca, 'FontSize', 36);

%% Linear regression computations: Traditional RSA estimates (Porges-Bohrer) vs. RSA amplitude from sinusoidal waves
x = [rsa_ibi_modelData.A]'; %(mean_acd_s0+mean_acd_s1+mean_acd_s2)./3;
% x = x(setdiff(1:23,find(isoutlier(baselineHFs_brSubj))));
y = [rsa_ibi_modelData.B]'; %(setdiff(1:23,find(isoutlier(baselineHFs_brSubj))));

x_c = x; %(~isoutlier(x) & ~isoutlier(y)); %there are no outliers in A, so I am discarding offset (B) outlier indices from both x and y.
y_c = y; %(~isoutlier(x) & ~isoutlier(y));

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
plot(x_c, y_c, 'o','MarkerSize',10,'Color','k','LineWidth',3);         % Plot original data points
hold on;
plot(x_c, y_fit, '-','LineWidth',3);     % Plot the fitted line
hold off;
xlabel('Mean ACD');
ylabel('B');
title('Mean ACD vs B');
legend('Data', 'Fitted Line');

% Add text to the northwest corner
text_loc_x = min(x_c); % Use the minimum x value for positioning
text_loc_y = max(y_c); % Use the maximum y value for positioning

% Calculate correlation coefficient and p-value
[r, p_value] = corr(x_c, y_c);

set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
text(text_loc_x, text_loc_y, ...
    sprintf('r = %.2f\np = %.4f', r, p_value), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline', ...
    'FontSize', 48, 'BackgroundColor', 'white', 'EdgeColor', 'black');
set(gca, 'FontSize', 36);


%% the amount of overall cardiac deceleration during stimulus anticipation (not taking breathing into account)
hit_dACD = mean_acd_s0_hit-mean_acd_s2_hit;
miss_dACD = mean_acd_s0_miss-mean_acd_s2_miss;

var1=hit_dACD;
% var1(11) = [];
var2=miss_dACD;
% var2(11) = [];

var12 = [{var1}; {var2}];

[p h] = signrank(var1,var2);

figure
h1 = rm_raincloud(var12,[0 0.8 0],0,'ks',[],0.01);

h1.p{2,1}.FaceColor=[0.8 0 0]; %face color for the second patch
h1.s{2,1}.MarkerFaceColor=[0.8 0 0]; %marker taskColor for the second category
h1.l(1,1).Visible="off";
h1.m(2,1).MarkerFaceColor = [0.8 0 0];

yticklabels({'Miss','Hit'}) %inverted because of the rm_raincloud function plot orientation!

% Also match individual data points
hold on
for i=1:size(var1,1)
    X = [h1.s{1,1}.XData(1,i) h1.s{2,1}.XData(1,i)];
    Y = [h1.s{1,1}.YData(1,i) h1.s{2,1}.YData(1,i)];
    plot(X,Y,"k-");
end

% Plot three lines as your significance line
% 1
hold on
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5];
Yv = [mean(h1.s{1,1}.YData) mean(h1.s{2,1}.YData)];
plot(Xv,Yv,"k-");
hold on
% 2
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+3 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5];
Yv = [mean(h1.s{1,1}.YData) mean(h1.s{1,1}.YData)];
plot(Xv,Yv,"k-");
% 3
Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+3 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+5];
Yv = [mean(h1.s{2,1}.YData) mean(h1.s{2,1}.YData)];
plot(Xv,Yv,"k-");

hold on
pF = round(p,2);
text(max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+10,(mean(h1.s{1,1}.YData)+mean(h1.s{2,1}.YData))/2,['p = ',num2str(pF)],"FontSize",36,"HorizontalAlignment","center",'FontWeight','bold');
xlim([-70 105])
xticks([-60 -40 -20 0 20 40 60 80 100])
set(gca,'FontSize',36,'FontWeight','Bold')
set(gcf, 'Position', get(0, 'Screensize'));
xlabel('Delta Offset') %we need to invert the x and y for all the labels and ticks as well due to rainplots being scripted in a rotated way!
title(['Total Cardiac Deceleration: Hit vs Miss'])
hold off

%% Diff between hit and miss
var1=hit_dACD;
% var1(11) = [];
var2=miss_dACD;
% var2(11) = [];
figure
raincloud_plot(var1 - var2,'box_on', 1,'color',[0.5 0.5 0.5]);
hold on
l1 = xline(0,'color',[0.2 0.2 0.2],'LineWidth',3,'LineStyle','--','Label','Origin','LabelHorizontalAlignment','left','LabelOrientation','horizontal','FontSize',32);
% l2 = xline(median(hit_avgPrePupil_mm - miss_avgPrePupil_mm),'color',[0.8 0 0],'LineWidth',3,'LineStyle','--');
yticks([]);
xlabel('Delta ACD Difference');
xlim([-40 45])
title(['Hit Minus Miss: Delta ACD'])
set(gca,'FontSize',36,'FontWeight','Bold')
xticks([-30 -20 -10 0 10 20 30 40])

%% plot the offset increase over the course of stimulus anticipation
for subj=1:23
    x2_group(subj,:) = rsa_ibi_modelData_BBB(subj).fitX_s2; y2_group(subj,:) = rsa_ibi_modelData_BBB(subj).fitY_s2;
    x1_group(subj,:) = rsa_ibi_modelData_BBB(subj).fitX_s1; y1_group(subj,:) = rsa_ibi_modelData_BBB(subj).fitY_s1;
    x0_group(subj,:) = rsa_ibi_modelData_BBB(subj).fitX_s0; y0_group(subj,:) = rsa_ibi_modelData_BBB(subj).fitY_s0;
    subjectNo = number_strings{1,subj};
    figure;
    plot(rsa_ibi_modelData_BBB(subj).fitX_s2, rsa_ibi_modelData_BBB(subj).fitY_s2, 'Color', [0.2 0.2 0.2], 'LineWidth', 3);
    hold on
    plot(rsa_ibi_modelData_BBB(subj).fitX_s1, rsa_ibi_modelData_BBB(subj).fitY_s1, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    hold on
    plot(rsa_ibi_modelData_BBB(subj).fitX_s0, rsa_ibi_modelData_BBB(subj).fitY_s0, 'Color', [0.8 0.8 0.8], 'LineWidth', 3);
    legend({'s-2','s-1','s-0'})
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 36)
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    figDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo];
    
    filename = [figDir '/offsetIncrease_s2_s1_s0_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [figDir '/offsetIncrease_s2_s1_s0_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;
end

%%% GROUP-LEVEL VERTICAL OFFSET PROGRESSION PLOT
figure;
plot(mean(x2_group,1),mean(y2_group,1),'Color', [0.2 0.2 0.2], 'LineWidth', 5);
hold on
plot(mean(x1_group,1),mean(y1_group,1),'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
plot(mean(x0_group,1),mean(y0_group,1),'Color', [0.8 0.8 0.8], 'LineWidth', 5);
legend({'s-2','s-1','s-0'})
xlabel('Breathing Phase (rad)')
ylabel('Delta IBI')
set(gca, 'FontSize', 36)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
figDir = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_figures_EK\GroupLevel\breathingDay1';
filename = [figDir '/offsetIncrease_s2_s1_s0_groupLevel.fig'];
saveas(gcf,filename)
filename = [figDir '/offsetIncrease_s2_s1_s0_groupLevel.png'];
saveas(gcf,filename)
close;
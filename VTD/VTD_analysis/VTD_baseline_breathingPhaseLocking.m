%% Do the IBI - Breathing phase Sinusoidal Fits for each individual
clear;clc
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\eeglab2023.0
eeglab

%% Get relevant file names from your subject of choice
subjectNo = input('Enter the code for the subject you want to analyze (e.g. 01): ', 's');
targetDir = strcat('S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_noneEEG_preProcessedBaselines\Cardiac_Breathing\',subjectNo);
filename = ['VTD' subjectNo '_ECG_BR_physEvents.set']; %timings of r-peaks, and breathing min (troughs) added!

EEG = pop_loadset('filename', filename, 'filepath', targetDir);
eeglab redraw
%% you need to obtain:
%%% 1) breathing phase of each R-peak
%%% 2) previous IBI and next IBI to calculate the delta IBI around that R-peak
srate = 1000;
card_and_br = EEG.event(strcmp({EEG.event.code},'b') | strcmp({EEG.event.code},'r'));
only_br = EEG.event(strcmp({EEG.event.code},'b'));
only_br_lats = [only_br.latency] .*(1000/srate);
only_card = EEG.event(strcmp({EEG.event.code},'r'));
only_card_lats = [only_card.latency] .*(1000/srate);

%% get the breathing phase of each R-peak
b=0;
for ev=1:size(only_card,2)
    b=b+1;
    rLat = only_card(ev).latency *(1000/srate);
    prevLats = only_br_lats(only_br_lats<rLat);
    postLats = only_br_lats(only_br_lats>=rLat);
    if isempty(prevLats) || isempty(postLats)
        brPhase(b,1) = NaN;
        brCycle(b,1) = NaN;
    else
        prevBR = max(prevLats);
        postBR = min(postLats);
        brCycle(b,1) = postBR - prevBR;
        brPhase(b,1) = ((rLat-prevBR)/(postBR-prevBR))*360; %assign a breathing phase for now, will check for breathing cycle duration criteria below
    end
    clear prevLats postLats
end
%% 
for br = 1:size(only_br,2)-1
    all_brCycles(br,1) = (only_br(br+1).latency - only_br(br).latency) *(1000/srate);
end
mean_brCycle = mean(all_brCycles);
std_brCycle = std(all_brCycles);

%%% Now that you know the mean and STD of breathing cycles for this subject, apply the exclusion criteria
for ph = 1:length(brPhase)
    if brCycle(ph,1) > mean_brCycle + (2 * std_brCycle)
        brPhase(ph,1) = NaN;
    elseif brCycle(ph,1) < mean_brCycle - (2 * std_brCycle)
        brPhase(ph,1) = NaN;
    end
end

%% get the previous and next IBI, to calculate the amount of ACD
for ev=1:size(only_card,2)
    rLat = only_card(ev).latency *(1000/srate);
    prevRs = only_card_lats(only_card_lats<rLat);
    postRs = only_card_lats(only_card_lats>rLat);
    if isempty(prevRs) || isempty(postRs)
        deltaIBI(ev,1) = NaN;
        deltaIBI(ev,1) = NaN;
    else
        prevIBI = rLat - max(prevRs);
        postIBI = min(postRs) - rLat;
        deltaIBI(ev,1) = postIBI - prevIBI;
    end
    clear prevRs postRs
end

%% now you have the brPhase and ACD (delta IBI) variables for each heartbeat.

%%% PERFORM THE RESTING STATE RSA CURVES (sinusoidal fit)
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_rsaSinusoidalModel
% Example data (replace with your actual vectors)
% phase_rad = ...     % Nx1 vector of breathing phase (radians)
% amplitude = ...     % Nx1 vector of response magnitude (e.g., ΔIBI)
brPhaseRad = deg2rad(brPhase);
brPhaseValid = brPhaseRad(~isnan(brPhase) & ~isnan(deltaIBI));
deltaIBIValid = deltaIBI(~isnan(brPhase) & ~isnan(deltaIBI));

[A, B, phi, y, SSres, SStot, Rsq, phase, amp, fitX, fitY] = RSA_sinusoidalFit(brPhaseValid,deltaIBIValid);

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_statisticsScripts\circstat-matlab-master\circstat-matlab-master
% Run circular-linear correlation
[rho, pval] = circ_corrcl(phase, amp);

% Display results
fprintf('Circular-linear correlation (actual delta IBI): rho = %.3f, p = %.4f\n', rho, pval);

% Visualization
figure;

% Plot as polar scatter
subplot(1,2,1)
polarplot(phase, amp, 'ko')
title('Polar Plot of Phase vs Delta IBI')
set(gca, 'FontSize', 14)

% Plot as linear scatter
subplot(1,2,2)
scatter(phase, amp, 40, 'k', 'filled')
hold on
plot(fitX,fitY, 'r-', 'LineWidth', 2)
%     hold on
%     plot(theta_fit_bias, amplitude_fit_bias, 'r--','LineWidth', 2)
xlabel('Breathing Phase (rad)')
ylabel('Delta IBI')
title([sprintf('Circular-linear Correlation\nrho = %.2f, p = %.4f', rho, pval) ' Subj: ' subjectNo]);
set(gca, 'FontSize', 14)

% Format the string
titleStr = sprintf('Real Data-- Fitted Amplitude: %.3f | Phase Shift: %.1f° | Offset: %.3f | R-sq: %.3f', ...
    A, rad2deg(phi), B, Rsq);

% Add it as the sgtitle
sgtitle(titleStr, 'FontSize', 16, 'FontWeight', 'bold');

%%% save the figure
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
filename = [targetDir '/circ2linear_corr_brPhase_ACD_baseline_subj' subjectNo '.fig'];
saveas(gcf,filename)
filename = [targetDir '/circ2linear_corr_brPhase_ACD_baseline_subj' subjectNo '.png'];
saveas(gcf,filename)
close;

%% you have the valid breathing phases and delta IBI:
%%% compute the average deltaIBI for 8 bins of breathing phase: [0-45], [45-90],...,[315 360] degrees (fog Figure 2B)
ibi_cat1 = []; ibi_cat2 = []; ibi_cat3 = []; ibi_cat4 = [];
ibi_cat5 = []; ibi_cat6 = []; ibi_cat7 = []; ibi_cat8 = [];

for p=1:length(deltaIBIValid)
    if brPhaseValid(p,1)>=deg2rad(0) && brPhaseValid(p,1)<deg2rad(45) %bin 1 [0 45]
        ibi_cat1 = [ibi_cat1 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(45) && brPhaseValid(p,1)<deg2rad(90) %bin 2
        ibi_cat2 = [ibi_cat2 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(90) && brPhaseValid(p,1)<deg2rad(135) %bin 3
        ibi_cat3 = [ibi_cat3 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(135) && brPhaseValid(p,1)<deg2rad(180) %bin 4
        ibi_cat4 = [ibi_cat4 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(180) && brPhaseValid(p,1)<deg2rad(225) %bin 5
        ibi_cat5 = [ibi_cat5 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(225) && brPhaseValid(p,1)<deg2rad(270) %bin 6
        ibi_cat6 = [ibi_cat6 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(270) && brPhaseValid(p,1)<deg2rad(315) %bin 7
        ibi_cat7 = [ibi_cat7 deltaIBIValid(p,1)];
    elseif brPhaseValid(p,1)>=deg2rad(315) && brPhaseValid(p,1)<deg2rad(360) %bin 8
        ibi_cat8 = [ibi_cat8 deltaIBIValid(p,1)];
    end
end
cat1_mean = mean(ibi_cat1,'omitnan'); cat2_mean = mean(ibi_cat2,'omitnan'); cat3_mean = mean(ibi_cat3,'omitnan'); cat4_mean = mean(ibi_cat4,'omitnan');
cat5_mean = mean(ibi_cat5,'omitnan'); cat6_mean = mean(ibi_cat6,'omitnan'); cat7_mean = mean(ibi_cat7,'omitnan'); cat8_mean = mean(ibi_cat8,'omitnan');

filename_rsa = 'rsa_effectBaseline_fig2B.mat';
    save(fullfile(targetDir, filename_rsa), "cat1_mean","cat2_mean","cat3_mean","cat4_mean","cat5_mean","cat6_mean","cat7_mean","cat8_mean");
%% compute the number of R-peaks in each breathing cycle: a metric for cardioventilatory coupling?
% Find points where the phase decreases (new cycle starts)
brPhase = brPhase(~isnan(brPhase));
cycleStartIdx = [1, (find(diff(brPhase) < 0) + 1)'];

% Append one beyond the last index to compute counts
cycleStartIdx(end+1) = length(brPhase) + 1;

% Count events per cycle
numCycles = length(cycleStartIdx) - 1;
eventsPerCycle = zeros(1, numCycles);

for i = 1:numCycles
    eventsPerCycle(i) = cycleStartIdx(i+1) - cycleStartIdx(i);
end

meanNum = mean(eventsPerCycle);
stdNum = std(eventsPerCycle);

for ev=1:size(eventsPerCycle,2)
    if (eventsPerCycle(i) > meanNum+3*stdNum) || (eventsPerCycle(i) < meanNum-3*stdNum)
        eventsPerCycle(i)=[];
    end
end
meanEvPerCycle= mean(eventsPerCycle); stdEvPerCycle = std(eventsPerCycle);

%% compute the traditional RSA estimation metrics:
%%% HF, pHF
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK\HRV
for c=2:size(only_card,2)
    ibi(c-1,1) = only_card(c).latency-only_card(c-1).latency;
end
ibiS = ibi ./1000;
[pLF,pHF,LFHFratio,VLF,LF,HF] = fft_val(ibiS,0,srate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Porges-Bohrer
% Step 1: Interpolate to get evenly sampled IBI
fs_interp = 4; % 4 Hz interpolation
t_ibi = cumsum(ibiS);
t_interp = t_ibi(1):1/fs_interp:t_ibi(end);
ibi_interp = interp1(t_ibi, ibi, t_interp, 'pchip');

% Step 2: Bandpass filter (typical adult RSA band: 0.12–0.40 Hz)
low_cutoff = 0.12; high_cutoff = 0.40;
[b,a] = butter(4, [low_cutoff high_cutoff] / (fs_interp/2), 'bandpass');
ibi_filtered = filtfilt(b, a, ibi_interp);

% Step 3: Sliding window RSA calculation (60-second window, 50% overlap)
win_length_sec = 60;
step_sec = win_length_sec;

win_samples = round(win_length_sec * fs_interp);
step_samples = round(step_sec * fs_interp);

n = length(ibi_filtered);
idx = 1:step_samples:(n - win_samples);

rsa_lnvar = zeros(length(idx),1);
for i = 1:length(idx)
    segment = ibi_filtered(idx(i):idx(i) + win_samples - 1);
    rsa_lnvar(i) = log(var(segment));
end

% Step 4: Plotting the RSA over time
time_rsa = t_interp(idx + floor(win_samples/2));
figure;
plot(time_rsa, rsa_lnvar, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('RSA (ln[ms²])');
title(['RSA via Porges-Bohrer Method: ' num2str(mean(rsa_lnvar)) 'ln[ms²]']);
grid on;
rsa_pb = mean(rsa_lnvar);
close;

disp(['RSA from sinusoidal: ' num2str(A)]);
disp(['RSA from HF: ' num2str(HF)]);
disp(['RSA from Porges-Bohrer: ' num2str(rsa_pb)]);
%% save the model params
filename_sineModel = 'rsa_ibi_sineModelBaseline.mat';
    save(fullfile(targetDir, filename_sineModel), 'A', 'phi','B','y','SSres','SStot','Rsq','HF','pHF','rsa_pb');
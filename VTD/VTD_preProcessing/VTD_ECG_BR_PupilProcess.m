% %% VTD_MAINProcess
% These set of scripts will be dealing with individual level pre-processing and...
% ...analyses
clear; clc;

%% Rename a file from Brainvision ONLY VIA THIS FUNCTION if necessary
% rename_brainvision_files('VTD28_Day1_Run3_Rest_restart.vhdr', 'VTD28_Day1_Run3_Rest.vhdr', 'rmf', 'on')
%% Get relevant file names from your subject of choice
subjectNo = input('Enter the code for the subject you want to analyze (e.g. 01): ', 's');

%% directs to target location where you have the ECG recordings of your...
%...subject
targetDir = strcat('S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_EEGrecordings\VTD',subjectNo);
cd(targetDir);

eegBlockFile = uigetfile({'*.*'},'Select the .eeg file of the block');
hdrBlockFile = uigetfile({'*.*'},'Select the .vhdr file of the block');

ecg1_chan = input('Enter the #ecg1 channel (most likely 65): ');
ecg2_chan = input('Enter the #ecg2 channel (most likely 66): ');

%% Start processing
% type in whether you want the following sanity checks:
% 1) Correct detection of R-peaks.
% 2) Correct detection of stimulus onset times within each trial.
visual = input('Do you want sanity check visualizations regarding correct phase and onset detections? (yes=1, no=0): ');

addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_analysisScripts
%% Read data in continuous mode!
eeg = char(eegBlockFile);
hdr = char(hdrBlockFile);
[cfgC, dataC] = continuousRead(eeg, hdr, subjectNo);
continuousRawECG1 = dataC.trial{1, 1}(ecg1_chan,:); % 1=ECG1 channel
continuousRawECG2 = dataC.trial{1, 1}(ecg2_chan,:); % 2=ECG2 channel

% if strcmp(subjectNo, '16') && contains(eegBlockFile, 'Day2_Run8', 'IgnoreCase', true)
%     subj16_weirdX = [277000,563200];
%     for el=1:length(subj16_weirdX)
%         continuousRawECG1(subj16_weirdX(el)-500 : subj16_weirdX(el)+500) = 0; %
%         dataC.trial{1,1}(1,subj16_weirdX(el)-500 : subj16_weirdX(el)+500) = 0;
%         dataC.trial{1,1}(2,subj16_weirdX(el)-500 : subj16_weirdX(el)+500) = 0;
%     end
%     continuousRawECG1(617300:619400) = 0; %
%     dataC.trial{1,1}(1,617300:619400) = 0;
%     dataC.trial{1,1}(2,617300:619400) = 0;
% end
%% detect R-peaks from the continuous mode and add the R-peak indices to the 'data' struct
if contains(eeg, 'Baseline' , 'IgnoreCase', true)
    [qrs_amp_raw, qrs_i_raw, delay, ecg_h] = pt_rpeaks(continuousRawECG1,1000,0); %Ash's script (last input argument is visualization: 1 for figures, 0 for no figures)

    for j=1:length(dataC.time{1,1})
        if ismember(j,qrs_i_raw)
            dataC.rPeaks{1,1}(1,j) = 1;
        else
            dataC.rPeaks{1,1}(1,j) = 0;
        end
    end
end
%% if you found the location of mistaken R-peak(s), correct it (them) here and restart the MAINProcess with the correction(s) uncommented!
dataC.rPeaks{1,1}(1,189736) = 0;
dataC.rPeaks{1,1}(1,188998) = 0;
% dataC.rPeaks{1,1}(1,1:4000) = 0;

%% sanity check plots for correct R-peak detection
%you can plot a random part of a trial a couple of times to see if...
%...the R-peak placement is correct throughout the continuous ECG.
%Choose 5 random numbers between [recordStart+10 recordEnd-30]
% (10 and 30 are seconds)

if contains(eeg, 'Baseline' , 'IgnoreCase', true) && visual == 1 %here we plot 5 random time intervals to sanity check that R-peaks are correctly detected.
    fs = dataC.fsample;
    N=5; %plot 5 random intervals
    r = 10*fs + ((length(dataC.time{1,1})-30*fs)-10000) .* rand(N, 1);
    plotLength = 10*fs; %10 second interval will be plotted.
    figure
    for f = 1:N
        subplot(N,2,f)
        plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),dataC.trial{1,1}(ecg1_chan,r(f):r(f)+plotLength).*-1); %raw ECG 

        hold on
        plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),dataC.rPeaks{1,:}(1,r(f):r(f)+plotLength) .* 1000); %time-points of R-peaks

        hold on
        plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),dataC.trial{1,1}(ecg2_chan,r(f):r(f)+plotLength));
%             hold on
%             plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),ecg_h(r(f)-1:r(f)+plotLength-1) .* 1000); %filtered (5-15 Hz bandpass for quality R-peak detection) ecg
        axis tight
        if f==1
            legend({'RawECG1','R-peaks','RawECG2'});
        end
        hold on
    end
    hold off
    %% Check the plausibility of interbeat intervals throughout the block (there should be NO outlier R-R interval value)
    rPeakPts = find(dataC.rPeaks{1,1});
    subplot(N,2,f+1)
    for i=1:length(rPeakPts)-1
        ibi(i) = rPeakPts(i+1) - rPeakPts(i);
        if ibi(i)>1500 %you can change this interval according to the subject, but not having a heartbeat in 2 seconds is not good :)
            ibi(i) = NaN;
%             dataC.trial{1,1}(rPeakPts(i):rPeakPts(i+1))=0;
        end
    end
    histogram(ibi/1000,20)
    sgtitle(['Sanity Check - R peak Detection ' eeg])
    rDetect_base = input('Are all R-peaks detected correctly (1=yes, 0=no)?: ');
end

if contains(eeg, 'Baseline' , 'IgnoreCase', true) && rDetect_base == 0
    plot(ibi); %to see where approximately there is a wrong IBI value!
    %% from the IBI plot, see which IBI is faulty and go to that area on the raw ECG plot:
    falseIBI = [211];
    for i=1:length(falseIBI)
        figure
        plot(dataC.time{1,1}(rPeakPts(falseIBI(i)-5):rPeakPts(falseIBI(i)+5)),dataC.rPeaks{1, 1}(1,(rPeakPts(falseIBI(i)-5):rPeakPts(falseIBI(i)+5))).*500);
        hold on
        plot(dataC.time{1,1}(rPeakPts(falseIBI(i)-5):rPeakPts(falseIBI(i)+5)),dataC.trial{1, 1}(ecg2_chan,rPeakPts(falseIBI(i)-5):rPeakPts(falseIBI(i)+5)).*-1);
%         plot(dataC.time{1,1},dataC.trial{1, 1}(1,:).*-1);
        hold off
    end
    return %stop running the code if there is an error in R-peak detection
elseif contains(eeg, 'Baseline' , 'IgnoreCase', true) && rDetect_base == 1
    dataC.ibiAll = ibi;
    close all
end


%% STIMULUS ONSET TIME POINT DETERMINATION
% Read trial structure from relevant files
if contains(eeg, 'Baseline' , 'IgnoreCase', true)
    %% get the average heart rate from the entire baseline recording
    dataC.avgIBI = mean(ibi,'omitnan');
    dataC.avgHeartRate = (1000/dataC.avgIBI)*60;

    addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_analysisScripts_EK\HRV 
    hrv_rmssd = RMSSD(ibi,0,0);
    hrv_sdsd = SDSD(ibi,0,0);
    hrv_sdnn = SDNN(ibi,0,0);
    hrv_pNN50 = pNN50(ibi,0,0);
    [pLF,pHF,LFHFratio,VLF,LF,HF] = fft_val(ibi,0,fs);
    hr = HR(ibi);
    hr = hr*1000; %because it was  /ms, now I am making it second.
    [med,qr,shift] = rrHRV(ibi);

    dataC.rmssd = hrv_rmssd;
    dataC.sdsd = hrv_sdsd;
    dataC.sdnn = hrv_sdnn;
    dataC.pNN50 = hrv_pNN50;
    dataC.pLF = pLF;
    dataC.pHF = pHF;
    dataC.LFHFratio = LFHFratio;
    dataC.VLF = VLF;
    dataC.LF = LF;
    dataC.HF = HF;
    dataC.avgHeartRate = hr;

    saveDirDesktop = strcat('C:\Users\kingir\Desktop\VTD_preProcessedBaselines\Cardiac_Breathing','\',subjectNo);
    mkdir(saveDirDesktop);
    cd(saveDirDesktop);
    blockNameSplit = split(eeg,'.');

    saveFileName = strcat(blockNameSplit{1}, '_cfgDataC');
    save(saveFileName, 'cfgC','dataC','-v7.3');

    saveDirDrive = strcat('S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_preProcessedBaselines\Cardiac_Breathing','\', subjectNo);
    mkdir(saveDirDrive);
    cd(saveDirDrive);
    save(saveFileName, 'cfgC','dataC','-v7.3');
    clear;clc;
    return %If the current block is a baseline, we will not go on with trial-based functions!
else
    [cfgTr, dataTr] = VTDtrialRead(eeg, hdr, subjectNo);
%     if size(dataC.label,1) == 68  
%         [dataTr] = photoDiodeFlashOnset(cfgTr, dataTr);
%     end
end

%% Determine the limits of trial period within the recording
firstTrialOnset = cfgTr.trl(1,1);
lastTrialEnd = cfgTr.trl(cfgTr.trialdef.nTR,2);
fs = dataC.fsample;

if contains(eeg, 'Rest' , 'IgnoreCase', true)
    %% Detect the R-peaks from continuous reading between the first trial onset and last trial's end
    [qrs_amp_raw, qrs_i_raw, delay, ecg_h] = pt_rpeaks(continuousRawECG1(firstTrialOnset:lastTrialEnd),1000,0); %Ash's script (last input argument is visualization: 1 for figures, 0 for no figures)
%     dataC.allRpeaks = ACCEPT_PHYSIO_MARK(ACCEPT_PHYSIO_MARK>=firstTrialOnset & ACCEPT_PHYSIO_MARK<lastTrialEnd);
    
    qrs_i_raw = qrs_i_raw + firstTrialOnset-1;
    for j=1:length(dataC.time{1,1})
        if ismember(j,qrs_i_raw)
            dataC.rPeaks{1,1}(1,j) = 1;
        else
            dataC.rPeaks{1,1}(1,j) = 0;
        end
    end

    rPeakPts = find(dataC.rPeaks{1,1});
    for i=1:length(rPeakPts)-1
        ibi(i) = rPeakPts(i+1) - rPeakPts(i);
        if ibi(i)>1500 %you can change this interval according to the subject, but not having a heartbeat in 2 seconds is not good :)
            ibi(i) = NaN;
%             dataC.trial{1,1}(rPeakPts(i):rPeakPts(i+1))=0;
        end
    end

    %% if you found the location of mistaken R-peak(s), correct it (them) here and restart the MAINProcess with the correction(s) uncommented!
%     dataC.rPeaks{1,1}(1,319072) = 0;
%     dataC.rPeaks{1,1}(1,242501) = 0;
%     dataC.rPeaks{1,1}(1,242582) = 1;
%     dataC.rPeaks{1,1}(1,242501) = 0;

    duringTrialsRpeaks = find(dataC.rPeaks{1,1}(firstTrialOnset:lastTrialEnd));
    for i=1:length(duringTrialsRpeaks)-1
        ibiDuringAllTrials(i) = duringTrialsRpeaks(i+1) - duringTrialsRpeaks(i);
    end
    dataTr.duringTrialsAvgHeartRate = (1000/mean(ibiDuringAllTrials,'omitnan'))*60;
    %% sanity check plots for correct R-peak detection
    %you can plot a random part of a trial a couple of times to see if...
    %...the R-peak placement is correct throughout the continuous ECG.
    %Choose 5 random numbers between [recordStart+10 recordEnd-30]
    % (10 and 30 are seconds)
    
    
    N=5; %plot 5 random intervals
    r = firstTrialOnset + 10*fs + ((length(firstTrialOnset:lastTrialEnd)-10000) .* rand(N, 1));
    plotLength = 10*fs; %10 second interval will be plotted.
    figure
    for f = 1:N
        subplot(N,2,f)
        plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),dataC.trial{1,1}(ecg1_chan,r(f):r(f)+plotLength).*-1); %raw ECG 

        hold on
        plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),dataC.rPeaks{1,:}(1,r(f):r(f)+plotLength) .* 1000); %time-points of R-peaks

        hold on
        plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),dataC.trial{1,1}(ecg2_chan,r(f):r(f)+plotLength));
%             hold on
%             plot(dataC.time{1,1}(1,r(f):r(f)+plotLength),ecg_h(r(f)-1:r(f)+plotLength-1) .* 1000); %filtered (5-15 Hz bandpass for quality R-peak detection) ecg
        axis tight
        if f==1
            legend({'RawECG1','R-peaks','RawECG2'});
        end
        hold on
    end
    hold off
    %% Check the plausibility of interbeat intervals throughout the block (there should be NO outlier R-R interval value)
    rPeakPts = find(dataC.rPeaks{1,1});
    subplot(N,2,f+1)
    for i=1:length(rPeakPts)-1
        if ibiDuringAllTrials(i)>1500 %you can change this interval according to the subject, but not having a heartbeat in 2 seconds is not good :)
            ibiDuringAllTrials(i) = NaN;
%             dataC.trial{1,1}(rPeakPts(i):rPeakPts(i+1))=0;
        end
    end
    dataC.allTrialsIBIs = ibiDuringAllTrials;
    histogram(ibiDuringAllTrials/1000,20)
    sgtitle(['Sanity Check - R peak Detection ' eeg])
    rDetect_rest = input('Are all R-peaks detected correctly (1=yes, 0=no)?: ');
    if rDetect_rest==0
        figure
        plot(ibiDuringAllTrials); %to see where approximately there is a wrong IBI value!
        hold off
        error('R-peaks not detected on a Rest block? What is wrong?');
    end
    

    %% from the IBI plot, see which IBI is faulty and go to that area on the raw ECG plot -- (default is commented)
    falseIBI = [141];
    for i=1:length(falseIBI)
        figure
        plot(dataC.time{1,1}(rPeakPts(falseIBI(i)-2):rPeakPts(falseIBI(i)+3)),dataC.rPeaks{1, 1}(1,(rPeakPts(falseIBI(i)-2):rPeakPts(falseIBI(i)+3))).*500);
        hold on
        plot(dataC.time{1,1}(rPeakPts(falseIBI(i)-2):rPeakPts(falseIBI(i)+3)),dataC.trial{1, 1}(ecg1_chan,rPeakPts(falseIBI(i)-2):rPeakPts(falseIBI(i)+3)).*-1);
%         plot(dataC.time{1,1},dataC.trial{1, 1}(1,:).*-1);
        hold off
    end

    %% Compute HRV params from the block (during trials)
    addpath C:\Users\kingir\Desktop\VTD_analysisScripts_EK\HRV 
    hrv_rmssd = RMSSD(ibi,0,0);
    hrv_sdsd = SDSD(ibi,0,0);
    hrv_sdnn = SDNN(ibi,0,0);
    hrv_pNN50 = pNN50(ibi,0,0);
    [pLF,pHF,LFHFratio,VLF,LF,HF] = fft_val(ibi,0,fs);
    hr = HR(ibi);
    hr = hr*1000; %because it was  /ms, now I am making it second.
    [med,qr,shift] = rrHRV(ibi);

    dataC.rmssd = hrv_rmssd;
    dataC.sdsd = hrv_sdsd;
    dataC.sdnn = hrv_sdnn;
    dataC.pNN50 = hrv_pNN50;
    dataC.pLF = pLF;
    dataC.pHF = pHF;
    dataC.LFHFratio = LFHFratio;
    dataC.VLF = VLF;
    dataC.LF = LF;
    dataC.HF = HF;
    dataC.avgHeartRate = hr;



    %% Detect relevant RR intervals within each trial
    for i =1:cfgTr.trialdef.nTR
        if contains(eegBlockFile, 'Rest', 'IgnoreCase', true)
            dataTr.trialBased_R_peaks{1,i} = dataC.rPeaks{1,1}(cfgTr.trl(i,1):cfgTr.trl(i,2));
        elseif contains(eegBlockFile, 'Cycle', 'IgnoreCase', true)
            dataTr.trialBased_R_peaks{1,i} = ACCEPT_PHYSIO_MARK(ACCEPT_PHYSIO_MARK>=cfgTr.trl(i,1) & ACCEPT_PHYSIO_MARK<cfgTr.trl(i,2));
        end
    end

    %% detect the phase of the stimulus within the cardiac and breathing cycles
    % cardiac phase
    %stimulus phase function determines the phase of stimulus onset within...
    %...cardiac cycle, and also assigns which bin the trial belongs to...
    %...when we categorize stimulus onset phases into bins:
    % 1) [0 200ms) wrt the previous R-peak
    % 2) [200 400ms)
    % 3) [400 600ms)
    % 4) [600 800ms)
    [dataTr] = VTDstimulusPhaseDetection(cfgTr, dataTr, cfgC, dataC); %creates a new struct in dataTr with stimulus onset phase values of each trial

    %% T-peak and offset detection
    % Necessary for "end of systole" detection, then systole and diastole...
    % ... assignment for eachY trial!
    tPeakDetected = 0;
    while tPeakDetected == 0
        [cfgTr, dataTr, invertData] = VTD_tPeakDetect(cfgC, dataC, cfgTr, dataTr,ecg1_chan,ecg2_chan);
        tPeakDetection = input('Are all t-peaks detected correctly? (1=yes, 0=no): ');
        if tPeakDetection == 1
            tPeakDetected = 1;
            close all
            break
        end
    end


    %% T-offset detection function also assigns each trial to systole, diastole, or in between.
    tOffsetDetected = 0;
    while tOffsetDetected == 0
        %Once t-peaks are detected, now you can go to the next step: Detecting...
        %...t-offsets
        [cfgTr, dataTr] = VTD_tOffsetDetect(cfgTr, dataTr, invertData, ecg1_chan, ecg2_chan, cfgC, dataC);
        tOffsetDetection = input('Are the detected t-offset locations correct (1=yes, 0=no)?: ');
        if tOffsetDetection == 1
            tOffsetDetected = 1;
            close all
            break
        end
    end

    %% Plot your target R-R intervals with filtered ECG, systole and diastole being shaded, and stimulus onset as an xline:
    figure
    t = tiledlayout(cfgTr.trialdef.nTR/6,cfgTr.trialdef.nTR/10,'TileSpacing','Compact');
    for i=1:cfgTr.trialdef.nTR
        nexttile
        if cfgTr.trl(i,3)==2
            continue
        
        elseif dataTr.tPeakErrors(i,1) == 0 && cfgTr.trl(i,3)~=2
            
            X1 = [dataTr.preRs(i) dataTr.preRs(i) find(dataTr.tOffsetPoint{1,i}) find(dataTr.tOffsetPoint{1,i})];
            Y1 = [0 500 500 0];
            p1 = patch(X1, Y1, 'r', 'FaceAlpha',.3);
            hold on
            X2 = [dataTr.postRs(i)-dataTr.systoleDiastoleLength(i) dataTr.postRs(i)-dataTr.systoleDiastoleLength(i) dataTr.postRs(i) dataTr.postRs(i)];
            Y2 = [0 500 500 0];
            p2 = patch(X2, Y2, 'b', 'FaceAlpha',.3);
            hold on

            stimulusOnset2Plot = cfgTr.trl(i,4)-cfgTr.trl(i,1);


            xline(stimulusOnset2Plot,'k--','LineWidth',1.5);
            targetTimeWindow = dataTr.preRs(i)-100:dataTr.postRs(i)+100;
            plot(targetTimeWindow,dataTr.trialBased_R_peaks{1, i}(targetTimeWindow).*500);
            hold on
            if dataTr.chan_tOffset==1
                plot(targetTimeWindow,dataTr.trialBased_tOffset_ecg1{1,i}(targetTimeWindow));
            elseif dataTr.chan_tOffset==2
                plot(targetTimeWindow,dataTr.trialBased_tOffset_ecg2{1,i}(targetTimeWindow).*-1);
            end
            xticks([]);
            axis tight
        end
    end
    
    xlabel(t,'Time limited by target IBI','FontSize',20);
    ylabel(t,'Raw ECG amplitude','FontSize',20);
    set(gcf, 'Position', get(0, 'Screensize'));
    hold off
    
    cardiacFigure = input('Are trial-based cardiac cycle characterizations good to save (1=yes, 0=no): ');

    c_fig_dir = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Cardiac'];
    if cardiacFigure==1
        blockNameSplit = split(eeg,'.');
        fig_dir = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Cardiac\'];
        mkdir(fig_dir);
        figure_png = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Cardiac\' blockNameSplit{1} '_cardiacTrials' '.png'];
        figure_mat = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Cardiac\' blockNameSplit{1} '_cardiacTrials' '.fig'];
        saveas(gcf,figure_png);
        saveas(gcf,figure_mat);
    end
    close all
    %% Define extended trials for breathing signal processing
    bb_chan = input('Enter the channel number for the breathing belt (most likely=67): ');
    for i =1:cfgTr.trialdef.nTR
        if i<cfgTr.trialdef.nTR
            dataTr.trialBased_breathing{1,i} = dataC.trial{1, 1}(bb_chan,(cfgTr.trl(i,1)-5*dataC.fsample):(cfgTr.trl(i,2)+5*dataC.fsample)); %adding 5 seconds to both ends of each trial, to plot the breathing activity better (because it is quite low freq.)
        elseif i==cfgTr.trialdef.nTR
            dataTr.trialBased_breathing{1,i} = dataC.trial{1, 1}(bb_chan,(cfgTr.trl(i,1)-5*dataC.fsample):(cfgTr.trl(i,2)+1*dataC.fsample));
        end
    end
    %% Detect IBI values in each trial: Will be used to perform similar analysis to Park et al. 2014 -- Figure 2: Post-decisional heart slowing
    [cfgTr, dataTr] = VTD_ibiDetect(cfgTr, dataTr,cfgC,dataC);
    
    %% Detect the phase of stimulus onset within the breathing cycle in each trial
    breathingPhaseDetected = 0;
    [cfgTr, dataTr, breathingSignalUsable] = VTD_breathingPhaseDetect(cfgTr, dataTr, dataC, visual);
    if breathingSignalUsable == 1
        breathingPhaseDetection = input('Are all breathing peaks detected correctly (1=yes, 0=no)?: ');
        if breathingPhaseDetection == 1
            breathingPhaseDetected = 1;
        end
    end


    if breathingSignalUsable==1 && breathingPhaseDetected==0
        invalidNo = input('Type the number of bad trials to be manually discarded: ');
        for i=1:invalidNo
            badTrial = input('Enter the next bad trial #: ');
            dataTr.stimBreathingPhases(1,badTrial) = NaN; %this means that this trial cannot be subjected to breathing phase analysis at later stages!
            dataTr.stimBreathingPhases_NoTrough(1,badTrial) = NaN;
            dataTr.stimBreathingQuadrant(1,badTrial) = NaN;
            dataTr.breathingCycleDuration(1,badTrial) = NaN;
        end

    end
    
    b_fig_dir = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Breathing'];
    mkdir(b_fig_dir);
    blockNameSplit = split(eeg,'.');
    figure_png = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Breathing\' blockNameSplit{1} '_breathingTrials' '.png'];
    figure_mat = ['C:\Users\kingir\Desktop\VTD_figures\' subjectNo '\Breathing\' blockNameSplit{1} '_breathingTrials' '.fig'];

    saveas(gcf,figure_png);
    saveas(gcf,figure_mat);
elseif contains(eegBlockFile, 'Cycl', 'IgnoreCase', true)
    addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\eeglab2023.0
    eeglab
    dataset_loaded = input('Type "1" when you successfully loaded the channel #65 from the relevant vhdr file to EEGLAB: ');

    EEG.data(1,:)=-1*EEG.data(1,:); %now the channel no is 1 because you upload only the ECG1 channel to the EEGLAB (write channel=65 on the pop-up menu)
    cardiac_events_accepted = input('Type "1" when you finish editing cardiac event markers on the EEGLAB GUI: ');
    dataC.allRpeaks = ACCEPT_PHYSIO_MARK(ACCEPT_PHYSIO_MARK>=firstTrialOnset & ACCEPT_PHYSIO_MARK<lastTrialEnd);
    
    for i=1:length(dataC.allRpeaks)-1
        ibiDuringAllTrials(i) = dataC.allRpeaks(i+1) - dataC.allRpeaks(i);
    end
    dataTr.duringTrialsAvgHeartRate = (1000/mean(ibiDuringAllTrials,'omitnan'))*60;

    %% Compute HRV params from the block (during trials)
    addpath C:\Users\kingir\Desktop\VTD_analysisScripts\HRV
    hrv_rmssd = RMSSD(ibiDuringAllTrials,0,0);
    hrv_sdsd = SDSD(ibiDuringAllTrials,0,0);
    hrv_sdnn = SDNN(ibiDuringAllTrials,0,0);
    hrv_pNN50 = pNN50(ibiDuringAllTrials,0,0);
    [pLF,pHF,LFHFratio,VLF,LF,HF] = fft_val(ibiDuringAllTrials,0,fs);
    hr = HR(ibiDuringAllTrials);
    hr = hr*1000; %because it was  /ms, now I am making it second.
    [med,qr,shift] = rrHRV(ibiDuringAllTrials,'omitnan');

    dataC.rmssd = hrv_rmssd;
    dataC.sdsd = hrv_sdsd;
    dataC.sdnn = hrv_sdnn;
    dataC.pNN50 = hrv_pNN50;
    dataC.pLF = pLF;
    dataC.pHF = pHF;
    dataC.LFHFratio = LFHFratio;
    dataC.VLF = VLF;
    dataC.LF = LF;
    dataC.HF = HF;
    dataC.avgHeartRate = hr;

    %% Detect relevant RR intervals within each trial
    for i=1:cfgTr.trialdef.nTR
        dataTr.trialBased_R_peaks{1,i} = dataC.allRpeaks(dataC.allRpeaks>=cfgTr.trl(i,1) & dataC.allRpeaks<cfgTr.trl(i,2));
        dataTr.tPeakErrors(i,1) = 1;
        dataTr.tPeak_tOffsetDistance(1,i) = NaN;
    end
    
    [dataTr] = VTDstimulusPhaseDetection(cfgTr, dataTr,cfgC,dataC);
    %% Detect IBI values in each trial: Will be used to perform similar analysis to Park et al. 2014 -- Figure 2: Post-decisional heart slowing
    [cfgTr, dataTr] = VTD_ibiDetect(cfgTr, dataTr,cfgC,dataC);
	
    %% Use this small section as a notebook to type in the mean systole lengths of each subject, from, the rest blocks of Day2
    % Subj. 01 = 308.5 --> 309
    % Subj. 02 = 334.3 --> 334
    % Subj. 03 = 339.0 --> 339
    % Subj. 04 = 324.5 --> 325
    % Subj. 05 = 316.6 --> 317
    % Subj. 06 = 328.1 --> 328
    % Subj. 07 = 298.6 --> 299
    % Subj. 08 = 333.7 --> 334
    % Subj. 09 = 295.1 --> 295
    % Subj. 10 = NaN because we are currently excluding the second day from this subject, as the eye-tracking fixation protocol did not work.
    % Subj. 11 = 344.8 --> 345
    % Subj. 12 = 357.3 --> 357
    % Subj. 13 = 345.4 --> 345
    % Subj. 14 = 337.6 --> 338
    % Subj. 15 = 282.1 --> 282
    % Subj. 16 = 328.3 --> 328
    % Subj. 17 = 349.5 --> 350
    % Subj. 18 = 341.8 --> 342
    % Subj. 19 = 309.5 --> 310
    % Subj. 20 = 303.3 --> 303
    % Subj. 21 = 347.8 --> 348
    % Subj. 22 = 312.4 --> 312
    % Subj. 23 = 308.1 --> 308
    % Subj. 24 = 313.9 --> 314
    % Subj. 25 = 331.0 --> 331
    % Subj. 26 = 370.4 --> 370
    % Subj. 27 = 314.1 --> 314
    % Subj. 28 = 335.9 --> 336

%     meanSys = mean(systoleMean, 'omitnan');
%     systoleMean = dataTr.systoleDiastoleLength;
%     systoleMean = [systoleMean dataTr.systoleDiastoleLength]


    systoleMean = input('Type in the closest integer to the mean systole length for this subject, as obtained from the Rest blocks on the same day (Day 2): ');
    for i=1:cfgTr.trialdef.nTR
        if cfgTr.trl(i,3) ~=2
            stimOnset = cfgTr.trl(i,4)-cfgTr.trl(i,1);
            systoleOnset = dataTr.preRs(1,i);
            systoleOffset = systoleOnset + systoleMean;
            diastoleOffset = dataTr.postRs(1,i);
            diastoleOnset = diastoleOffset - systoleMean;
            if stimOnset>=systoleOnset && stimOnset<=systoleOffset
                dataTr.systole_or_diastole(1,i) = 1; %1=systole
            elseif stimOnset>=diastoleOnset && stimOnset<=diastoleOffset
                dataTr.systole_or_diastole(1,i) = 2; %2=diastole
            else
                dataTr.systole_or_diastole(1,i) = 0; %0=neither systole nor diastole
            end
        end
    end
    
    

end

%% categorize Hit vs. Miss from the MATLAB recordings
dayNo = input('Type whether the recording is from day 1 or day 2 of the subject (type 1 or 2): ');
runNo = input('Type the block number (between 1 and 8): ');

behaviorDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_MATLABRecordings\VTD' subjectNo '\Day' num2str(dayNo)];
cd(behaviorDir);
behFile = ['VTD_Subj' subjectNo '_Run' num2str(runNo)];

behList = dir(fullfile(behaviorDir));
behFileNames = behList(~cellfun(@(x) any(regexp(x, '^\.+$')), {behList.name})); % avoid '.' and '..'

for beh = 1:length(behFileNames)
    if strfind(behFileNames(beh).name,behFile)
        blockBehavior = load(behFileNames(beh).name);
    end
end
hitOrMiss = NaN(cfgTr.trialdef.nTR,1);

for tr=1:cfgTr.trialdef.nTR
    if blockBehavior.output(tr+10).condition == 1 && blockBehavior.output(tr+10).response == 1 %we put +10 because we are discarding the first 10 trials from each block.
        hitOrMiss(tr,1) = 1;
    elseif blockBehavior.output(tr+10).condition == 1 && blockBehavior.output(tr+10).response == 0
        hitOrMiss(tr,1) = 0;
    elseif blockBehavior.output(tr+10).condition == 2
        hitOrMiss(tr,1) = -1; %-1=stimulus absent trial
    end
end

%% make a good matrix structure that contains all relevant information for the initial merged-block analysis step
time = dataTr.time';
trial = dataTr.trial';
rPeaks = dataTr.trialBased_R_peaks';
stimulusPhase = dataTr.stimulusPhase';
stimulusBin200 = dataTr.stimulusBin200';
stimulusBin50 = dataTr.stimulusBin50';
stimulusQuadrant = dataTr.stimulusQuadrant';
hitMiss = hitOrMiss;
trialIBIs = dataTr.trial_IBIs';
day = repmat(dayNo,cfgTr.trialdef.nTR,1);
run = repmat(runNo,cfgTr.trialdef.nTR,1);
type    = cell(cfgTr.trialdef.nTR,1);
type(:) = {blockBehavior.params.runtype};

for i=1:cfgTr.trialdef.nTR
    contrastMean(i) = blockBehavior.output(i+10).AllContrastEstimatesMean; %hard-coded between 11:70 for now, means that I am skipping the first 10 trials in a 70-trial block.
    contrastQuantile(i) = blockBehavior.output(i+10).AllContrastEstimatesQuantile;
    confidence(i) = blockBehavior.output(i+10).confidenceAnswer;
end
contrastMean = contrastMean';
contrastQuantile = contrastQuantile';
confidence = confidence';

if contains(eegBlockFile, 'Rest', 'IgnoreCase', true)
    systoleOrDiastole = dataTr.systole_or_diastole_tOffsetMethod;
    breathingPhase = dataTr.stimBreathingPhases';
    breathingPhaseNoTrough = dataTr.stimBreathingPhases_NoTrough';
    breathingQuadrant = dataTr.stimBreathingQuadrant';
    blockTable = table(day,run,type,time,trial,rPeaks,stimulusPhase,stimulusBin200,stimulusBin50,stimulusQuadrant,hitMiss,contrastMean,contrastQuantile,confidence,systoleOrDiastole,trialIBIs,breathingPhase,breathingPhaseNoTrough,breathingQuadrant);

elseif contains(eegBlockFile, 'Cycl', 'IgnoreCase', true)
    systoleOrDiastole = dataTr.systole_or_diastole';
    blockTable = table(day,run,type,time,trial,rPeaks,stimulusPhase,stimulusBin200,stimulusBin50,stimulusQuadrant,hitMiss,contrastMean,contrastQuantile,confidence,systoleOrDiastole,trialIBIs);
end
    

%% save the cfgC, dataC, cfgTr, dataTr structs to a convenient location
% C=continuous (R-peaks)
% Tr=trial-based (stimulus phase and T-wave analyses etc.)

blockNameSplit = split(eeg,'.');
saveFileName = strcat(blockNameSplit{1}, '_cfgDataStructs');
saveDirDesktop = strcat('C:\Users\kingir\Desktop\VTD_preProcessed\Cardiac_Breathing','\',subjectNo); %Save to the S-drive
mkdir(saveDirDesktop);
cd(saveDirDesktop);
save(saveFileName, 'cfgC','dataC','cfgTr','dataTr','blockTable','-v7.3');

saveDirDrive = strcat('S:\KNEU\KNEUR-Projects\Projects\Sukanya-Backup\VTD_preProcessed\Cardiac_Breathing','\', subjectNo);
mkdir(saveDirDrive);
cd(saveDirDrive);
save(saveFileName, 'cfgC','dataC','cfgTr','dataTr','blockTable','-v7.3');

clear;clc



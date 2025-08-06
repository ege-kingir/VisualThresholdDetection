%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VTD - EEG Pre-Processing %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test for ASV file function

%% Overall Preprocessing Scheme

% ================================================================================================================= %
% IN THIS PROJECT, EEG ANALYSIS IS PERFORMED AFTER ECG-BB-PUPIL ANALYSES SO WE ALREADY HAVE THE TRIAL STRUCTURE...
% ================================================================================================================= %

% ================================================================================================================= %
% STRONG RECOMMENDATION:
% To avoid errors, and also keep track of your pre-processing steps; run this code Section-by-Section!
% ================================================================================================================= %

% Pre-processing Steps:

% 1a) Open EEGLAB, load the EEG data (only channels 1:64), do this for each block before anything else; ...
% ... and have all raw datasets corresponding to each block from the subject loaded in EEGLAB.

% 2a) Load ECG-BB data from each block and get the R-Peaks and determined Stimulus Onset time points (1 block from a subject).
% ... sections 2a and 2b should be run in order for each block ...
% ... you need to switch the active dataset on EEGLAB into whichever block events you are obtaining ...

% 2b) BEFORE FILTERING: Concatenate all blocks temporarily to see which channels look deviant in the power spectra...
% ... (plot spectrum between 2-80 hz). Note the "absolutely" deviant channels for removal LATER!
% Also visualize via channel scroll to verify that these channels are visibly bad on time-domain?

% 3a) Low-pass and high-pass filter through Basic FIR filter in EEGLAB (bandpass between 0.1 and 40 hz)...
% 3b) CleanLine for line noise removal: only change to default settings is the first one ("Freq to remove"): make it [50 100] instead of 60
% ... Save continuous block sets at 1000 hz (no downsampling yet): "vtdXX_d1rX_1000hzCont"

% 4) Go back to separated block datasets (Make sure you saved the continuous data above, you are gonna use those again for downsampling).
% ...EPOCH between ([0 6.001] with respect to 'S  1' and 'S  2':
%       Gets all valid trials in a block, and 6.001 instead of 6 is for Fieldtrip compatibility.)
% ...MUSCLE artifact rejection through Fieldtrip...
% ...Run section 4a for each block.
% ...Section 4 in the script also gets info from PupilLabs data to detect trials that include blinks (via PupilLabs Core analysis)...
% ...See the variable "isblink" in your workspace: Take note of blink-trials because you have to reject them. They are blinks that happened in...
% ...our critical time window ([-1.5 +0.3] seconds with respect to stimulus onset.

% ...Most z-scored artifacts are actually benign (Fieldtrip magnifies the effect of high-freq. activity: You can select "really bad epochs" from these).
% ...Take note of trials that deviate "a lot (trials that show a visible level of noise in many electrodes)" via manual judgment. But do not reject them yet.

% ...You have to do the step 4 above because you can use Fieldtrip only with 1000 Hz data (in accordance with your cfg structure)...
% ...But you want to proceed with 500 Hz data, so you perform step 5a and onwards as described below:

% 5a) Go back to EEGLAB GUI and load continuous block sets for each block again...
% ... Downsample the continuous data to 500 Hz...
% ... Epoch all of the blocks with respect to S1 and S2, between [0 6.001] sec
% ... Then discard the artifactual/blink trials (the ones you took note of).

% ... EPOCH 1: Concatenate the blocks and convert to stimulus-aligned blocks in order to better visualize "stimulus-related" events on ICA components:
%       Epoch into [-3 1.5] seconds with respect to 'S  3', which is the stimulus onset event trigger
%       There might be a few trials that give this error on the CMD: "Warning: event XXX out of data boundary"...
%       You need to add these trials to the rejected ones because there was an error in the trigger timing for sitmulus onset.
%       Once done, save as "...preICA_s3..."

% 5b) Run section 5b on this script to get the MATLAB behavioral values (hit-miss, breathing phase, cardiac phase, eye-tracking values...
% ... and confidence levels etc.)
% ... SAVE the workspace variable "eegValidTrialStruct" when you get valid trials from each block.
% ... Also check that the total # of trials (epochs) are equal on "eegValidTrialStruct" (sum the # of rows taken for each block) and EEGLAB.
% ... NOTE: Valid trials should be exactly the same for both EPOCH 1 and EPOCH 2

% ... EPOCH 2: Convert the concatenated epochs to trial-start aligned again, but just locked to 'S 1' [0 6.001], discarding the stimulus absent ('S 2') trials.
%     Once done, save as "...preICA_s1..."

% 6a) Load two copies of the preICA_s3 dataset. Filter the second one with 1-30hz
% 6b) Run ICA (extended runica) but exclude the bad channel(s) that are bad in all or most blocks on the ICA menu.

% 7) Visualize components first by "inspect/label components by map" and then via ICLabel.
% 8) Look at the first 35 components and definitely exclude components that are >90% non-brain (but >90% a single category) according to ICLabel:
%       If you already noted a lot of frontal channels with muscle noise, you can remove the components that definitely look like muscle noise...
%       ... even when the ICLabel does not classify it as >90% single other category. These bad components should have no prominent alpha power...
%       ... and should have visibly increased beta power, hinting towards quite significant muscle noise. And the components should be...
%       ... mostly localized to frontal areas because the most prominent muscle noise source is straining the eyebrows. Some subjects do it a lot!

% 9) Decide on the components to be removed, note them down... 
%    ...Copy the ICA weights (no removal yet) to the light-filtered data (section 6a and 6b in the script).
%    ...Save both the original EEG (light filters) and the heavily filtered copy for ICA (you should make sure that both have ICA weights).
%       Then remove the components from the original EEG with ICA weights.
%    ...SAVE the "original EEG with bad components removed" as well.

% GO BACK TO STEP 6a, and this time repeat until this point with the S1-merged dataset...
%    Save your "badCompRem" dataset for the S1-merged one...
%    Also save the post ICA version of the data filtered for ICA.
%    ...Then do the necessary steps for cross-correlation with ECG signal to remove ECG components.

% THE STEPS BELOW SHOULD BE DONE FOR BOTH S1 and S3 MERGED DATASETS:

% 10) Interpolate if there are any bad channels that were excluded from ICA.
% IMPORTANT: Also save the "original EEG with bad comps removed + bad chans interpolated"

% 11) Go to section 8 of the script, load relevant datasets on EEGLAB and run the subsections accordingly to...
% ... save the valid trials that are only Hits, Misses, Low confidence, High confidence...
% ... two different subsections on the script are for separate categorization of visual ERP and HEP analysis versions of the datasets!

% NOTE: You can also re-reference to the common average and save that version as well (or you can do it from the dataset from step 12...
% ...any time you want).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               FOR HEARTBEAT EVOKED POTENTIAL ANALYSIS               %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% UPDATE 01.10.2024: I am switching to Pairwise Phase Consistency method...
% ...for a more reliable, literature-based detection of ECG-artifactual components.

%%% STEPS TO FOLLOW %%%
%%% with the EEG set:
% Take the dataset with bad components removed, but no interpolation of the bad channels yet! (s1Merged_filt4ica_badCompRem)

% Cross-Correlation: Epoch that wrt the R-peaks: [-0.2 0.35] seconds
% Pairwise Phase Consistency: Epoch that wrt the R-peaks: [-0.1 0.35] seconds
% Save the epoched dataset as "vtdXX_filt4ica_eeg_rPeakEpochs"

%%% with the ECG set:
% Load channel 65 data from each block (OR channel 66: look for the level of noise, decide which ECG channel represents a typical ECG better).
% Run section 2a of the script for each block to extract relevant parts of the recording, and with all valid trial events.

% Low- and high-pass the ECG data between 1 and 30 hz (changed from 0.1-40 to 1-30 because I decided to use the heavy filtered data at [1 30] for EEG as well).
% Apply CleanLine to each block.
% Downsample to 500hz and SAVE continuous data from each block.

% Epoch all blocks according to 'S 1' and 'S 2', then reject same epochs as you did for the EEG pipeline above.

% Append the datasets and epoch again, with respect to 'S 1': [0 6.001] seconds

% Yet another epoch again, now wrt the R-peaks (check that you obtain same amount of epochs with the EEG set epoched around R-peaks).
% The epochs should be [-0.2 0.35] if you are doing cross-correlation, [-0.1 0.35] if you are doing PPC.

% Save the ECG set epoched around the r-peaks: "vtdXX_rPeaksPPC_ecg" (the PPC denotes the switch to PPC method for ECG artifact detection)

% Run the "Cross-Correlation between ICA and ECG..." section of this code. This will give you the "ECG components to be removed"

% Go back to your EEG dataset: "original EEG with bad components removed", and remove the ECG components determined via Cross-Corr.:
    % IMPORTANT: Maximum number of components you can remove = 3
    % If you see more than 3 ECG comps, look at the amounts of correlations and decide on the top 3 high-priority components to remove.

% You can look at the effects of removing the ECG components, using BrainVision Analyzer visualizations:
    % Export the two datasets "original EEG with bad comps removed + bad chans interpolated"...
    % ... and "ECG comps removed + bad chans interpolated" separately to the Analyzer.
    % Then re-epoch both sets wrt the R-peaks within the "segmentation" tool in Analyzer.
    % Take averages (HEP) of both, and drag one on top of the other on the visualization window.


%% UPDATE:

%%%% The PREP step below, might be harmful for ICA because it references the data to the "estimate" of the true AVERAGE, but there may still be some
% noisy channels, which then may contaminate all other channels. So this step is taken out of pre-processing for now (02.05.24)

% 4) Through the EEGLAB GUI, run the PREP-pipeline which consists of...:
%       1) High-pass filter at 1 Hz. (THIS IS TEMPORARY, SO YOU WILL HIGH-PASS YOUR DATA AGAIN FROM 0.5 Hz AFTER PREP).
%       2) Remove line noise: Change the line noise frequencies from 60 to 50 and its multiples.
%       3) Referencing to the estimated "true" average EEG signal
%       4) Interpolating the bad channels due to: Extreme amplitudes, lack of correlation w/other channels, and too much high freq. noise:
%           BUT we are not interpolating now, we choose the box "remove bad interpolated channels" from the Post-process box of the PREP GUI!
%           Reason for no interpolation right now: We want the interpolation to be done after the data is cleaned via ICA, so that we can use...
%           "clean" data to interpolate the bad channel!
%           %%% It is worth noting that this early interpolation step does not remove physiological artifacts!

%%%% Instead of cross-correlation method to detect ECG-like artifactual ICA components, I switched to Pairwise Phase Consistency.


%% 1) Open EEGLAB
clear; clc
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\eeglab2023.0
eeglab

eegValidTrialStruct = []; %will be filled in section 5b


%%% add Fieldtrip path as well, will be necessary in the future
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\fieldtrip-20230613

folder = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD_EEGrecordings';
addpath(genpath(folder));

%% 2a) Load the cfgTr.trl from the ECG-BB analysis of the corresponding block...

subjectNo = input('Enter the code for the subject you want to analyze (e.g. 01): ', 's'); %VTD subject for analysis
day = input('Which day is the recording from (1 or 2): ');
runNo = input('Which block number is it (between 1 and 8): ');

targetDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_preProcessed\Cardiac_Breathing\' num2str(subjectNo)];
fileList = dir(fullfile(targetDir, ['*Day' num2str(day) '_Run' num2str(runNo) '_*'])); %fullfile=directs full file names from parts
blockStruct = load(fullfile(targetDir, fileList.name));

trialOrg = blockStruct.cfgTr.trl; %the matrix with start, end, stimulus onset time points of the trial!

allStimOnsets = trialOrg(:,4);
hitStimOnsets = trialOrg(trialOrg(:,3)==1,4); %Stimulus Onset Time Points from the trials where the subject successfully detected the faint stimulus.
missStimOnsets = trialOrg(trialOrg(:,3)==0,4); %Stimulus Onset Time Points from the trials where the subject missed the faint stimulus.

rPeaks = find(blockStruct.dataC.rPeaks{1, 1})'; %Time points of each R-peak, as detected previously during ECG-BB analysis!

% 2b) Filter the events such that you only have the ones that are important for you!
loaded = input('Type 1 on the CMD once you downloaded the dataset: '); %this way we wait until the dataset is loaded, to play around with the events on MATLAB...

%%% INDICES of the events that are in your valid trial structure: %%%
idx = ismember([EEG.event.latency],trialOrg(:,1)) | ismember([EEG.event.latency],trialOrg(:,4)) | ismember([EEG.event.latency],trialOrg(:,5));
idxNo = find(idx);
EEG.event = EEG.event(idxNo);

%%% ADD R-Peak events %%%
for i=1:length(rPeaks)
    EEG.event(end+1) = EEG.event(length(EEG.event)); %create a new event, by copying the last valid event to the next row as well.
    EEG.event(end).latency = rPeaks(i);
    EEG.event(end).type = 'r peak';
    EEG.event(end).code = 'cardiac event';
end

firstTr_start = round(trialOrg(1,1)/1000,3)
lastTr_end = round(trialOrg(end,2)/1000,3)

EEG = pop_select(EEG,'time',[firstTr_start-1 lastTr_end+1]);

if strcmp(subjectNo, '03')
    EEG.chanlocs = misplaceChanLocs; % I HAD TO DO THIS FOR SUBJECT 03 BECAUSE TWO SPLITTER BOXES ARE CONFUSED DURING DATA COLLECTION IN THAT CASE !!!
end


[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % Store dataset once you are done removing & adding events!

clearvars subjectNo day runNo targetDir fileList trialOrg allStimOnsets hitStimOnsets missStimOnsets rPeaks firstTr_start lastTr_end idx idxNo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CURRENTLY I CANNOT ADD EXPIRATION ONSETS FROM BREATHING DATA, BECAUSE I HAVE NOT SAVED THEM...%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4a) Muscle artifact detection via Fieldtrip (decisions will be made manually)

%%%% Epoch your EEG data into trials: epoch based on the 'S  1' and 'S  2' events, and window is [0 8] with respect to each 'S  1' or 'S  2' event.
    %%% NOTE: You do not do baseline subtraction for now!
    %%% Convert the EEG struct to a fieldtrip struct!
% EEG_ft_epoch = eeglab2fieldtrip(EEG,'raw','none'); %convert to fieldtrip structure after epoching the trials within EEGLAB (60 trials should be converted typically).

%%%% define the trials for fieldtrip (we will feed in the outputs that we had from ECG-BB analysis)

% ====================================================================================== %
            %%%% YOU NEED TO CONFIGURE THE CFG AND YOUR DATA PROPERLY %%%%
% Need to make sure that the structs are appropriate with necessary fields for fieldtrip %
% ====================================================================================== %



%%% Load the current block's ECG-BB data for trial structure extraction before Fieldtrip artifact detection

subjectNo = input('Enter the code for the subject you want to analyze (e.g. 01): ', 's'); %VTD subject for analysis
day = input('Which day is the recording from (1 or 2): ');
runNo = input('Which block number is it (between 1 and 8): ');

targetDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_preProcessed\Cardiac_Breathing\' num2str(subjectNo)];
fileList = dir(fullfile(targetDir, ['*Day' num2str(day) '_Run' num2str(runNo) '_*'])); %fullfile=directs full file names from parts
blockStruct = load(fullfile(targetDir, fileList.name));

%%%% You want to reject trials with blinks in the pre-stimulus and/or stimulus response window, so look at "isblink" to choose trials to reject.
expTr = find(blockStruct.blockTable.hitMiss ~= -1);
blink_and_catch = find(isnan(blockStruct.blockTable.avgPreStim_pupilSizeNorm) | isnan(blockStruct.blockTable.avgStimResponse_pupilSizeNorm));
isblink = intersect(expTr, blink_and_catch);

clearvars artifact_muscle artifact_jump cfg

%%% Prepare the structure for Fieldtrip compatibility %%%

%%% Adapting the original data to the pre-processing steps that I did within EEGLAB:
blockStruct.dataTr.trial = [];
if blockStruct.cfgTr.trialdef.trialLength==8
    blockStruct.dataTr.sampleinfo(:,2) = blockStruct.dataTr.sampleinfo(:,2)-2000;
end

for trial =1:size(EEG.data,3)
    blockStruct.dataTr.trial{1,trial} = EEG.data(:,:,trial); %decreases the data into the "accepted" channels from the pre-processing in EEGLAB until now (we are still pre-ICA).
end
blockStruct.dataTr.label = {EEG.chanlocs.labels};

for tr=1:size(blockStruct.dataTr.time,2)
    blockStruct.dataTr.time{1,tr} = blockStruct.dataTr.time{1,tr}(1:6001);
end


folder = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\rawData\VTD_EEGrecordings\VTD' subjectNo];
cd(folder);
addpath(genpath(folder));

channelExc = input('Type in the excluded channel numbers before this step: ');
allChans = 1:64;
channelInc = ~ismember(allChans,channelExc);

channelSel = blockStruct.dataC.label(channelInc);

%%% Muscle Artifacts via Fieldtrip z-score method %%%

% ============================================================================================================================================= %
% Look here to see default params for other artifact detections:
% https://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection
% ============================================================================================================================================= %

cfg = [];
cfg.trialdef.ntrials = size(EEG.data,3);
cfg.trl = blockStruct.cfgTr.trl(:,1:2);

if blockStruct.cfgTr.trl(1,2) - blockStruct.cfgTr.trl(1,1) ==8000
    cfg.trl(:,2) = cfg.trl(:,2)-2000;
end

cfg.sampleinfo = blockStruct.dataTr.sampleinfo;
cfg.continuous = 'no';
cfg.dataset = blockStruct.cfgC.dataset;
cfg.trialdef.triallength = 6;

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel      = channelSel;
cfg.artfctdef.zvalue.cutoff       = 4; %default was 4
cfg.artfctdef.zvalue.trlpadding   = -0.5; %negative trial padding is to just save the edges of trials from filter artifacts...
cfg.artfctdef.zvalue.fltpadding   = 0.5;
cfg.artfctdef.zvalue.artpadding   = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter     = 'yes';
cfg.artfctdef.zvalue.bpfreq       = [110 140];
% cfg.artfctdef.zvalue.bpfiltord    = 9;
cfg.artfctdef.zvalue.bpfilttype   = 'firws';
cfg.artfctdef.zvalue.hilbert      = 'yes';
cfg.artfctdef.zvalue.boxcar       = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';
% 
[cfg, artifact_muscle] = ft_artifact_zvalue(cfg,blockStruct.dataTr);

%%%% After looking at the possible artifact segments/channels for muscle artifacts, look into your data from ft_databrowser!

cfg.trl(:,3) = cfg.trl(:,2) - cfg.trl(:,1);
[cfg] = ft_databrowser(cfg,blockStruct.dataTr);

%% 5) For each block: List the rejected trials and get the MATLAB recording info for only accepted trials:
%%% eegValidTrialStruct will become important for categorizing EEG epochs according to low or high confidence, different breathing phases etc etc...

subjectNo = input('Enter the code for the subject you want to analyze (e.g. 01): ', 's'); %VTD subject for analysis
day = input('Which day is the recording from (1 or 2): ');
runNo = input('Which block number is it (between 1 and 8): ');

targetDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_preProcessed\Cardiac_Breathing\' num2str(subjectNo)];
fileList = dir(fullfile(targetDir, ['*Day' num2str(day) '_Run' num2str(runNo) '_*'])); %fullfile=directs full file names from parts
blockStruct = load(fullfile(targetDir, fileList.name));

trialOrg = blockStruct.cfgTr.trl; %the matrix with start, end, stimulus onset time points of the trial!

nTrials = size(trialOrg,1);
rejectedTrials = input('List rejected trials (in []): ');
stimAbsTrials = find(blockStruct.blockTable.hitMiss==-1)';
stimPresent2Keep = [rejectedTrials stimAbsTrials];
allTrials = [1:nTrials];
trials2keep = setdiff(allTrials,stimPresent2Keep);
eegValidTrialStruct{runNo} = blockStruct.blockTable(trials2keep,:);

%%% check the total number of stimulusPresent and Valid trials:
trs = 0;
for b=1:size(eegValidTrialStruct,2)
    trs = trs+size(eegValidTrialStruct{1,b},1);
end

%% EXTRA: Get the latest stimulus onset following the trial start via S1
% presents=0;
% for tr=1:EEG.trials
%     stimPresent = find(strcmpi(EEG.epoch(tr).eventtype,{'S  1'}));
%     if ~isempty(stimPresent)
%         presents=presents+1;
%         stimOnset = find(strcmpi(EEG.epoch(tr).eventtype,{'S  3'}));
%         stimLatency(presents) = EEG.epoch(tr).eventlatency{1, stimOnset};
%     end
% end
% 
% maxLatency = max(stimLatency);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   POST ICA                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6a) name 2 EEG files properly (AFTER YOU RUN ICA):

lightEEG = ALLEEG(1); %the one that we will do the actual analysis (ICA weights from the other set will be pasted here)
heavyEEG = ALLEEG(2); %the one which has the ICA weights (heavily filtered for this purpose only)

%%% get ica weights to temp variables
% out_icaact = heavyEEG.icaact;
out_icawinv = heavyEEG.icawinv;
out_icaweights = heavyEEG.icaweights;
out_icasphere = heavyEEG.icasphere;
out_icachansind = heavyEEG.icachansind;

%%% put these weights to the temp variable for the light-filtered EEG set
% lightEEG.icaact = out_icaact;
lightEEG.icawinv = out_icawinv;
lightEEG.icaweights = out_icaweights;
lightEEG.icasphere = out_icasphere;
lightEEG.icachansind = out_icachansind;

%% 6b) Here, make sure that the active EEG set on EEGLAB is the one with light filters (the one with no ICA weights yet).
EEG = lightEEG;

% EEG.chanlocs = misplaceChanLocs; % I HAD TO DO THIS FOR SUBJECT 03 BECAUSE TWO SPLITTER BOXES ARE CONFUSED DURING DATA COLLECTION IN THAT CASE !!!

[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % Store dataset once you are done removing & adding events!


%% 7a) 14.10.24 -- Adding ECG as an additional channel before ICA decomposition (with the filt4ica_badCompRem dataset)
% EEG_set = ALLEEG(1); %make sure that your dataset #1 is the EEG dataset, and epoched in [-0.1 0.35] around the R-peak
% ECG_set = ALLEEG(2); %this is the ECG dataset, epoched in the same way (and pre-processed before hand in the same manner as the EEG set)
% 
% EEG_dat = EEG_set.data;
% ECG_dat = ECG_set.data;
% 
% EEG_dat(end+1,:,:) = ECG_dat;
% EEG_set.data = EEG_dat;
% EEG_set.nbchan = EEG_set.nbchan + 1;
% EEG_set.chanlocs(end+1).labels = 'ecg';  % Label the new channel
% 
% EEG_set = eeg_checkset(EEG_set);  % Check the consistency of the dataset
% ALLEEG(1) = EEG_set;
% eeglab redraw;  % Update the EEGLAB GUI



%% 7b) Pairwise Phase Consistency between ECG and ICA-component activations to discard ICA components that show max. correlation values of >mean+3*SD with the ECG data.

% clear;clc;
% eeglab

EEG_set = ALLEEG(1); %make sure that your dataset #1 is the EEG dataset, and epoched in [-0.1 0.35] around the R-peak
ECG_set = ALLEEG(2); %this is the ECG dataset, epoched in the same way (and pre-processed before hand in the same manner as the EEG set)

EEG_dat = EEG_set.icaact;
ECG_dat = ECG_set.data;

filter=1; %for 
[ecg_artComps, PPC_allComps, PPC_wMaxLag, ppc_limit] = ppc_vtd_with_lags(EEG_dat, ECG_dat, filter);

%% 8 ) Categorize trials according to Hit-Miss and Low or High confidence:

%%% WARNING:
% For visual ERP purposes: Make sure that the active EEG dataset on EEGLAB is "vtdXX_badCompRem_interp"

%%  8a) For visual ERP analysis
clear; clc
addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\eeglab2023.0
eeglab

% Before running next section, For visual ERP purposes: Make sure that the active EEG dataset on EEGLAB is "vtdXX_badCompRem_interp"

%%
subjectNo = input('Enter the code for the subject you want to analyze (e.g. 01): ', 's');
validTrLoc = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\rawData\VTD_EEG_preprocessing\vtd' subjectNo];
cd(validTrLoc);
load(['vtd' subjectNo '_d1_stimPresentValidTrials']);

allBlocks = [];
for sess=1:size(eegValidTrialStruct,2)
    if ~isempty(eegValidTrialStruct{1,sess})
        allBlocks = [allBlocks; eegValidTrialStruct{1,sess}(:,1:19)];
    end
end

%%% a) Hit-Miss

hits = allBlocks(allBlocks.hitMiss==1,:);
hitInd = find(allBlocks.hitMiss==1);

miss = allBlocks(allBlocks.hitMiss==0,:);
missInd = find(allBlocks.hitMiss==0);

%%% b) Low-High Confidence (Median split)

medianConf = median(allBlocks.confidence);
lowConf = allBlocks(allBlocks.confidence<medianConf,:);
lowConfInd = find(allBlocks.confidence<medianConf);

highConf = allBlocks(allBlocks.confidence>medianConf,:);
highConfInd = find(allBlocks.confidence>medianConf);

hitEEG = pop_select(EEG,'trial',hitInd);
ALLEEG(end+1) = hitEEG;

missEEG = pop_select(EEG,'trial',missInd);
ALLEEG(end+1) = missEEG;

highConfEEG = pop_select(EEG,'trial',highConfInd);
ALLEEG(end+1) = highConfEEG;

lowConfEEG = pop_select(EEG,'trial',lowConfInd);
ALLEEG(end+1) = lowConfEEG;

savePath = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing\vtd' subjectNo];
pop_saveset(hitEEG,['vtd' subjectNo '_hitsVisERP'],savePath);
pop_saveset(missEEG,['vtd' subjectNo '_missesVisERP'],savePath);
pop_saveset(lowConfEEG,['vtd' subjectNo '_lowConfVisERP'],savePath);
pop_saveset(highConfEEG,['vtd' subjectNo '_highConfVisERP'],savePath);

%% 8b) For HEP analysis (temporary automatization for group-level extraction: 07.10.24)
clear;clc
conf_analysis = input('Are you analysing smth related to confidence levels (1:yes, 0:no): ');
breathing_analysis = input('Are you analysing smth related to breathing (1:yes, 0:no): ');

% Define the directory path containing the folders
folder_path = 'S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing_newFilters'; % Replace with your directory path

% Get a list of all folders in the directory
folder_list = dir(fullfile(folder_path, 'vtd*')); % Assuming folders start with 'vtd'

% Initialize an empty cell array to store the numeric parts as strings
number_strings = {};
% Loop through the folder list to extract the numeric parts
for i = 1:length(folder_list)
    folder_name = folder_list(i).name; % Get the folder name
    % Use regular expression to extract the numeric part
    num_str = regexp(folder_name, '\d+', 'match');
    if ~isempty(num_str)
        % Store the numeric part as a string with leading zeros (if necessary)
        formatted_num_str = sprintf('%02d', str2double(num_str{1}));
        number_strings = [number_strings, formatted_num_str];
    end
end

nSubj=27;
if conf_analysis==1
    nSubj=23; %subj 06 is excluded due to too few valid EEG trials in day1, and 4 more subjects are excluded from confidence analysis
    conf_delete = cellfun(@(x) contains(x, {'05','16','26','28'}), number_strings);
    number_strings(conf_delete) = [];
end

if breathing_analysis==1
    nSubj=23; %subj 06 is excluded due to too few valid EEG trials in day1, and 4 more subjects are excluded from confidence analysis
    conf_delete = cellfun(@(x) contains(x, {'03','08','20','24'}), number_strings);
    number_strings(conf_delete) = [];
end

conf_table50 = table('Size',[nSubj 2],'VariableTypes',{'double','double'},'VariableNames',{'below50','above50'});
highConf_hitRate = zeros(nSubj,1);
lowConf_hitRate = zeros(nSubj,1);
highConf_stimContrast = zeros(nSubj,1);
lowConf_stimContrast = zeros(nSubj,1);
hitConf = zeros(nSubj,1);
missConf = zeros(nSubj,1);
hitStimContrast = zeros(nSubj,1);
missStimContrast = zeros(nSubj,1);

for subj=1:nSubj
%     clear; clc
    subjectNo = number_strings{1,subj};
    addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\eeglab2023.0
    eeglab
%     mkdir(['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_epochs4LIMO\vtd' subjectNo]);
    savePath = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing_newFilters\vtd' subjectNo];
    limoPath = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_epochs4LIMO\vtd' subjectNo];
    loadPath = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing_newFilters\vtd' subjectNo];
    directory = fullfile(loadPath);
    fileListD = dir(fullfile(directory, ['*newFilt_s3_visualERP*.set'])); %04.02.25 --> changing to s3epochs because I want to use Visual ERP modeling with LIMO toolbox.
%     fileListD = dir(fullfile(directory, ['*s1epochs_point5_40_ecgCompRem*reref*.set']));
    [EEG] = pop_loadset(fileListD.name,loadPath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % Store dataset once you are done removing & adding events!
    
    % Before running next section, For HEP analysis purposes: Make sure that the active EEG dataset on EEGLAB is "vtdXX_ecgCompRem_interp"
    %%% UPDATE 09.09.24: I decided that I will first try without removing the <3 ECG-like components, so now I am running the section...
    %%% ...below with ECG components included.
    % So the dataset to use: "vtdXX_s1_badCompRem_interp" (if interpolation was not necessary, then use "vtdXX_s1_badCompRem").
    % And the section below saves the categories of Hit, Miss, LowConf, HighConf as "vtdXX_*category*HEP_s1_wECG"
    
    %%%

    
    validTrLoc = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing\vtd' subjectNo];
    cd(validTrLoc);
    load(['vtd' subjectNo '_d1_stimPresentValidTrials']);
    
    
    allBlocks = [];
    for sess=1:size(eegValidTrialStruct,2)
        if ~isempty(eegValidTrialStruct{1,sess})
            allBlocks = [allBlocks; eegValidTrialStruct{1,sess}(:,1:19)];
        end
    end    
    
    

%     if breathing_analysis==1
%         lateExhales = allBlocks(allBlocks.breathingPhaseNoTrough>90 & allBlocks.breathingPhaseNoTrough<180,:);
%         lateExhaleInd = find(allBlocks.breathingPhaseNoTrough>90 & allBlocks.breathingPhaseNoTrough<180);
%         
%         lateExhaleEEG = pop_select(EEG,'trial',lateExhaleInd);
%         ALLEEG(end+1) = lateExhaleEEG;
%         disp(['Late exhale trial count: ' num2str(length(lateExhaleInd))]);
%     end
    
    inhales=[];
    exhales=[];
    for i=1:length(EEG.epoch)
        EEG.epoch(i).cardPhase = cell2mat([allBlocks.stimulusPhase(i)]);

        if allBlocks.stimulusQuadrant(i) == 1
            EEG.epoch(i).cardQ = 'one';
        elseif allBlocks.stimulusQuadrant(i) == 2
            EEG.epoch(i).cardQ = 'two';
        elseif allBlocks.stimulusQuadrant(i) == 3
            EEG.epoch(i).cardQ = 'three';
        elseif allBlocks.stimulusQuadrant(i) == 4
            EEG.epoch(i).cardQ = 'four';
        end
        
        EEG.epoch(i).brPhase = allBlocks.breathingPhaseNoTrough(i);

        if allBlocks.breathingQuadrant(i)==1 || allBlocks.breathingQuadrant(i)==2
            EEG.epoch(i).InhaleOrExhale = 'exhale';
            exhales = [exhales i];
        elseif allBlocks.breathingQuadrant(i)==3 || allBlocks.breathingQuadrant(i)==4
            EEG.epoch(i).InhaleOrExhale = 'inhale';
            inhales = [inhales i];
        end
        
        EEG.epoch(i).detection = logical(allBlocks.hitMiss(i));
        EEG.epoch(i).confidence = allBlocks.confidence(i);
        EEG.epoch(i).contrast = allBlocks.contrastMean(i);
    end

    meanSurr_vep_str = load([savePath '\vtd' subjectNo '_meanSurr_VEP_fromMinus1toPlus1_combined']); %loading the mean surrogate VEP data
    meanSurr_vep = meanSurr_vep_str.meanSurr_vep;  
    
    EEG_s3 = pop_epoch(EEG,'S  3',[-1 1.002],'newname','stim_epochs');
    for tr=1:size(EEG_s3.data,3) %this is subtracting the mean surrogate VEP from the non-corrected VEP
        EEG_s3.data(:,:,tr) =  EEG_s3.data(:,:,tr) - meanSurr_vep;
        EEG_s3.epoch(tr).cardPhase = EEG.epoch(tr).cardPhase;
        EEG_s3.epoch(tr).cardQ = EEG.epoch(tr).cardQ;
        EEG_s3.epoch(tr).brPhase = EEG.epoch(tr).brPhase;
        if ~ismember(subjectNo,{'03','08','20','24'})
            EEG_s3.epoch(tr).InhaleOrExhale = EEG.epoch(tr).InhaleOrExhale;
        end
        EEG_s3.epoch(tr).detection = EEG.epoch(tr).detection;
        EEG_s3.epoch(tr).confidence = EEG.epoch(tr).confidence;
        EEG_s3.epoch(tr).contrast = EEG.epoch(tr).contrast;
    end

%     if ~ismember(subjectNo,{'03','08','20','24'})
%         exhaleEEG = pop_select(EEG,'trial',exhales);
%         inhaleEEG = pop_select(EEG,'trial',inhales);
%     
%         pop_saveset(exhaleEEG,['vtd' subjectNo '_newFilt_HEP_s1_exhaleTrials'],savePath);
%         pop_saveset(inhaleEEG,['vtd' subjectNo '_newFilt_HEP_s1_inhaleTrials'],savePath);
%     end
    
    %%% a) Hit-Miss
%     if conf_analysis==0 && breathing_analysis==0
%         hits = allBlocks(allBlocks.hitMiss==1,:);
%         hitInd = find(allBlocks.hitMiss==1);
%         
% % 
% %         [hit_minCont, min_ind] = sortrows(hits,'contrastQuantile','ascend'); %min_ind sorts the indices within the hit trials according to stimulus contrast
% % 
%         miss = allBlocks(allBlocks.hitMiss==0,:);
%         missInd = find(allBlocks.hitMiss==0);
% %         
% %         hit_minInds = min_ind(1:numel(missInd)); %downsamples the hit trials to the size of miss trials, and includes the hit trials with smallest possible stimulus contrast values while doing that.
% %         hit_minInds = sort(hit_minInds,'ascend');
% %         hitMinTable = hit_minCont(1:100,:);
% 
% %         hitConf(subj,1) = mean(hits.confidence);
% %         missConf(subj,1) = mean(miss.confidence);
% % 
%         hitStimContrast(subj,1) = mean(hits.contrastMean);
%         missStimContrast(subj,1) = mean(miss.contrastMean);
% 
% %         trialIBIs = cell2mat(allBlocks.trialIBIs);
% %         decelInd = find(trialIBIs(:,1)<trialIBIs(:,2) & trialIBIs(:,2)<trialIBIs(:,3));
% 
% %         hitEEG_ds = pop_select(hitEEG,'trial',hit_minInds);
%         
% %         decelEEG = pop_select(EEG,'trial',decelInd);
% 
%         hitEEG = pop_select(EEG,'trial',hitInd); %has all of the hit trials
%         ALLEEG(end+1) = hitEEG;
% % %         
%         missEEG = pop_select(EEG,'trial',missInd);
%         ALLEEG(end+1) = missEEG;
%     
%     %%% b) Low-High Confidence
%     elseif conf_analysis==1 && breathing_analysis==0
%         medianConf = median(allBlocks.confidence);
%         lowConf = allBlocks(allBlocks.confidence<=medianConf,:);
%         lowConfInd = find(allBlocks.confidence<=medianConf);
%         
%         highConf = allBlocks(allBlocks.confidence>medianConf,:);
%         highConfInd = find(allBlocks.confidence>medianConf);
%         
% %         if size(highConfInd,1)>size(lowConfInd,1)
% %             highConfDS = randsample(highConfInd,size(lowConfInd,1)); %take a random sample from high-confidence trials, according to the number of low confidence trials...
% %             highConfDS = sort(highConfDS,'ascend');
% %             highConfEEG = pop_select(EEG,'trial',highConfDS);
% %             lowConfEEG = pop_select(EEG,'trial',lowConfInd);
% %         elseif size(highConfInd,1)<size(lowConfInd,1)
% %             lowConfDS = randsample(lowConfInd,size(highConfInd,1));
% %             lowConfDS = sort(lowConfDS,'ascend');
% %             lowConfEEG = pop_select(EEG,'trial',lowConfDS);
% %             highConfEEG = pop_select(EEG,'trial',highConfInd);
% %         end
% 
%         lowConfEEG = pop_select(EEG,'trial',lowConfInd);
%         highConfEEG = pop_select(EEG,'trial',highConfInd);
%         highConf_hitRate(subj,1) = sum(highConf.hitMiss) / size(highConf,1);
%         lowConf_hitRate(subj,1) = sum(lowConf.hitMiss) / size(lowConf,1);
% 
%         highConf_stimContrast(subj,1) = mean(highConf.contrastMean);
%         lowConf_stimContrast(subj,1) = mean(lowConf.contrastMean);
% 
%         ALLEEG(end+1) = lowConfEEG;
%         ALLEEG(end+1) = highConfEEG;
%     end
    
%     if breathing_analysis==1
%         pop_saveset(lateExhaleEEG,['vtd' subjectNo '_newFilt_lateExhaleEEG_HEP_s1_noECG_reref'],savePath); %ds means that the trial numbers are downsampled to have equal trials in hit and miss conditions!
% %         pop_saveset(missEEG,['vtd' subjectNo '_newFilt_missesHEP_s1_noECG_reref'],savePath);
%     end
        
%     if conf_analysis==0 && breathing_analysis==0
%         pop_saveset(decelEEG,['vtd' subjectNo '_newFilt_decelEEG_HEP_s1_noECG_reref'],savePath);
% %         pop_saveset(hitEEG,['vtd' subjectNo '_newFilt_hits_NoDownSampling_HEP_s1_noECG_reref'],savePath); %ds means that the trial numbers are downsampled to have equal trials in hit and miss conditions!
% %         pop_saveset(missEEG,['vtd' subjectNo '_newFilt_missesHEP_s1_noECG_reref'],savePath);
%     elseif conf_analysis==1 && breathing_analysis==0
%         % LR = Left-Right = high confidence trials are trials where the subject slides the confidence bar to the Right, and vice versa but counting...
%         % ...  mid-point as well for low confidence trials.
%         pop_saveset(lowConfEEG,['vtd' subjectNo '_newFilt_medianSplit_lowConfHEP_s1_noECG_reref'],savePath); % LR = high confidence trials are trials where the subject slides the confidence bar to the right, and vice versa but counting mid-point as well for low confidence trials
%         pop_saveset(highConfEEG,['vtd' subjectNo '_newFilt_medianSplit_highConfHEP_s1_noECG_reref'],savePath);
%     end
%     pop_saveset(EEG,['vtd' subjectNo '_newFilt_visERP_combinedEpochs_forLIMO'],limoPath);
    pop_saveset(EEG_s3,['vtd' subjectNo '_newFilt_visERP_m1p1_100surrCorr_combinedEpochs_forLIMO'],limoPath); %this means that the VEP curves are also corrected for the non-stim-locked ongoing activity (just as i did for HEPs): m1p1 means the epochs start from -1 and go until +1 second wrt the target stimulus onset!
%     pop_saveset(EEG,['vtd' subjectNo '_newFilt_HEP_s1_combinedEpochs_forLIMO'],savePath);
    
    clearvars -except number_strings conf_analysis breathing_analysis highConf_hitRate lowConf_hitRate hitConf missConf highConf_stimContrast lowConf_stimContrast hitStimContrast missStimContrast
end

%% Compare mean stimulus contrasts in High Conf. vs. Low Conf. Trials
% var1=highConf_stimContrast;
% var2=lowConf_stimContrast;
% 
% var12 = [{var1}; {var2}];
% 
% p_stim = signrank(var1,var2);

%% Compare mean stimulus contrast in Hit (downsampled to include Hit trials with the lowest stimulus contrast values) vs. miss trials
% var1=hitStimContrast;
% var2=missStimContrast;
% 
% var12 = [{var1}; {var2}];
% 
% p_stimHitMiss = signrank(var1,var2);
% 
%% Plot Hit Rates in Low vs. High Confidence Trials
% addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\ViolinPlot\raincloudPlots
% var1=highConf_hitRate;
% var2=lowConf_hitRate;
% 
% var12 = [{var1}; {var2}];
% 
% p = signrank(var1,var2);
% 
% figure
% 
% h1 = rm_raincloud(var12,[0 0 0.8],0,'ks',[],1);
% 
% h1.p{2,1}.FaceColor=[0 0.8 0]; %face color for the second patch
% h1.s{2,1}.MarkerFaceColor=[0 0.8 0]; %marker taskColor for the second category
% h1.l(1,1).Visible="off";
% h1.m(2,1).MarkerFaceColor = [0 0.8 0];
% 
% yticklabels({'Low Confidence','High Confidence'}) %inverted because of the rm_raincloud function plot orientation!
% 
% % Also match individual data points
% hold on
% for i=1:size(var1,1)
%     X = [h1.s{1,1}.XData(1,i) h1.s{2,1}.XData(1,i)];
%     Y = [h1.s{1,1}.YData(1,i) h1.s{2,1}.YData(1,i)];
%     plot(X,Y,"k-");
% end
% 
% % Plot three lines as your significance line
% % 1
% hold on
% Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2];
% Yv = [mean(h1.s{1,1}.YData) mean(h1.s{2,1}.YData)];
% plot(Xv,Yv,"k-");
% hold on
% % 2
% Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2];
% Yv = [mean(h1.s{1,1}.YData) mean(h1.s{1,1}.YData)];
% plot(Xv,Yv,"k-");
% % 3
% Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2];
% Yv = [mean(h1.s{2,1}.YData) mean(h1.s{2,1}.YData)];
% plot(Xv,Yv,"k-");
% 
% hold on
% pF = round(p,3);
% text(max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.25,(mean(h1.s{1,1}.YData)+mean(h1.s{2,1}.YData))/2,['p = ',num2str(pF)],"FontSize",36,"HorizontalAlignment","center",'FontWeight','bold');
% xlim([0 1.2])
% xticks([0.2 0.4 0.6 0.8 1])
% set(gca,'FontSize',36,'FontWeight','Bold')
% set(gcf, 'Position', get(0, 'Screensize'));
% xlabel('Hit Rates') %we need to invert the x and y for all the labels and ticks as well due to rainplots being scripted in a rotated way!
% title(['Hit Rates in Low and High Confidence Trials'])
% hold off
% 
% %% Plot Confidence in Hits (downsampled) and Misses
% 
% addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\ViolinPlot\raincloudPlots
% var1=hitConf/100;
% var2=missConf/100;
% 
% var12 = [{var1}; {var2}];
% 
% p = signrank(var1,var2);
% 
% figure
% 
% h1 = rm_raincloud(var12,[0 0 0.8],0,'ks',[],1);
% 
% h1.p{2,1}.FaceColor=[0 0.8 0]; %face color for the second patch
% h1.s{2,1}.MarkerFaceColor=[0 0.8 0]; %marker taskColor for the second category
% h1.l(1,1).Visible="off";
% h1.m(2,1).MarkerFaceColor = [0 0.8 0];
% 
% yticklabels({'Misses','Hits'}) %inverted because of the rm_raincloud function plot orientation!
% 
% % Also match individual data points
% hold on
% for i=1:size(var1,1)
%     X = [h1.s{1,1}.XData(1,i) h1.s{2,1}.XData(1,i)];
%     Y = [h1.s{1,1}.YData(1,i) h1.s{2,1}.YData(1,i)];
%     plot(X,Y,"k-");
% end
% 
% % Plot three lines as your significance line
% % 1
% hold on
% Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2];
% Yv = [mean(h1.s{1,1}.YData) mean(h1.s{2,1}.YData)];
% plot(Xv,Yv,"k-");
% hold on
% % 2
% Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2];
% Yv = [mean(h1.s{1,1}.YData) mean(h1.s{1,1}.YData)];
% plot(Xv,Yv,"k-");
% % 3
% Xv = [max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.1 max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.2];
% Yv = [mean(h1.s{2,1}.YData) mean(h1.s{2,1}.YData)];
% plot(Xv,Yv,"k-");
% 
% hold on
% pF = round(p,3);
% text(max(max(h1.s{1,1}.XData,h1.s{2,1}.XData))+0.25,(mean(h1.s{1,1}.YData)+mean(h1.s{2,1}.YData))/2,['p = ',num2str(pF)],"FontSize",36,"HorizontalAlignment","center",'FontWeight','bold');
% xlim([0 1.3])
% xticks([0.2 0.4 0.6 0.8 1])
% set(gca,'FontSize',36,'FontWeight','Bold')
% set(gcf, 'Position', get(0, 'Screensize'));
% xlabel('Confidence') %we need to invert the x and y for all the labels and ticks as well due to rainplots being scripted in a rotated way!
% title(['Confidence in Hits (downsampled) and Misses'])
% hold off

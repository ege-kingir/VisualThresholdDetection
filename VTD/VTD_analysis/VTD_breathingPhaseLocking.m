%% VTD -- breathing phase locking analysis to each R-peak

% 07.02.25 --> I categorize Inhalation vs Exhalation based on breathing phase at stimulus onset
% (instead of the previous approach where categorize according to the breathing phase at each heartbeat separately)
% So in this new approach, an "exhalation" trial can have a breathing phase corresponding to inhalation at S-2, or S-1.

%%
clear; clc
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

nSubj=23; %subj 06 is excluded due to too few valid EEG trials in day1, and 4 more subjects are excluded from confidence analysis
conf_delete = cellfun(@(x) contains(x, {'03','08','20','24'}), number_strings);
number_strings(conf_delete) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NEEDS TO BE RUN INDIVIDUAL BY INDIVIDUAL !!! %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You unfortunately need to adjust event filters for some subjects individually, if the heart rate is higher in a subject, you need to adjust which 'S  3' indices are eligible for being the true stimulus indices for that subject...
    %%% ... you need to do the above mentioned adjustment only for breathing analysis, because for breathing peak detection you need extended trials,
    %%% ... which are computationally a bit troublesome.
    %%% Examples:
    %%% VTD07 (subj #5) --> if i==260; stimInd=7;
    %%% VTD10 (subj #7) --> ev<11
    %%% VTD15 (subj #12) --> ev<10 and uncomment this in the script:
%         if i==361
%             stimInd =  8;
%             epochEvents{i,stimInd} = 'stim'; epochEvents{i,10} = 'WrongStim'; epochEvents{i,23} = 'WrongStim';
%         end
    %%% VTD23 (subj #19) --> ev>=8 and uncomment this in the script
%         if i==310
%             stimInd =  7;
%             epochEvents{i,stimInd} = 'stim'; epochEvents{i,21} = 'WrongStim';
%         elseif i==208
%             stimInd = 8;
%             epochEvents{i,stimInd} = 'stim';
%         elseif i==335
%             epochEvents{i,9} = 'WrongStim';
%         end
    %%% VTD25 (subj=20) --> ev>=10 && ev<34 and uncomment this in the script
%                 if i==109
%                     currHit = 0;
%                     respLatency = epochLatency{i,34};
%                     respTime = respLatency - epochLatency{109,24};
%                     missRespTimes = [missRespTimes respTime];
%                     allRespTimes = [allRespTimes respTime];
%                 elseif i==348
%                     currHit = 1;
%                     respLatency = epochLatency{i,38};
%                     respTime = respLatency - epochLatency{348,29};
%                     hitRespTimes = [hitRespTimes respTime];
%                     allRespTimes = [allRespTimes respTime];
%                 end
for subj=1:1
    %% load the fast and slow response trials
    subjectNo = number_strings{1,subj};
    resps = load(['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo '\fast_slowResponseTrials.mat']);

    %% load the subject table that includes the confidence values too!
    behavTableLoc = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_preProcessed\CombinedCardiacBreathing\VTD' subjectNo '\1'];
    subjTable = load([behavTableLoc '\VTD' subjectNo '_Day1_combined_IBIplus3Incl']);

    if ~ismember(subjectNo,{'05','16','26','28'})
        confTable = subjTable.combined_withIBIplus3(subjTable.combined_withIBIplus3.hitMiss<2,:);
        medianConf = median(confTable.confidence,'omitnan');
        confs = load(['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo '\lowHighConfTrials.mat']);
        highC_trs = confs.highConfTrials;
        lowC_trs = find(confs.lowConfTrials);
    else
        highC_trs = [];
        lowC_trs = [];
    end
    %% NOTE from 14.02.2025: When you add troughs on top of peaks in breathing, the event numbers should change accordingly:
    % For example on average it seems like 4 breathing events are added. I might change the default ev>7 and ev<=24 to ev>11 and ev<=28

    exhale_s2 = [];
    exhaleHit_s2 = [];
    exhaleMiss_s2 = [];
    inhale_s2 = [];
    inhaleHit_s2 = [];
    inhaleMiss_s2 = [];
    exhale_s1 = [];
    exhaleHit_s1 = [];
    exhaleMiss_s1 = [];
    inhale_s1 = [];
    inhaleHit_s1 = [];
    inhaleMiss_s1 = [];
    exhale_s0 = [];
    exhaleHit_s0 = [];
    exhaleMiss_s0 = [];
    inhale_s0 = [];
    inhaleHit_s0 = [];
    inhaleMiss_s0 = [];
    exhaleHit_stim = [];
    exhaleMiss_stim = [];
    inhaleHit_stim = [];
    inhaleMiss_stim = [];
    exhale_s_p1 = [];
    exhaleHit_s_p1 = [];
    exhaleMiss_s_p1 = [];
    inhale_s_p1 = [];
    inhaleHit_s_p1 = [];
    inhaleMiss_s_p1 = [];
% 
%     accel_s2 = []; accelHit_s2=[]; accelMiss_s2=[]; accel_s2_phaseShift=[]; accel_s2_avgIBI=[];
%     decel_s2 = []; decelHit_s2=[]; decelMiss_s2=[]; decel_s2_phaseShift=[]; decel_s2_avgIBI=[];
%     accel_s1 = []; accelHit_s1=[]; accelMiss_s1=[]; accel_s1_phaseShift=[]; accel_s1_avgIBI=[];
%     decel_s1 = []; decelHit_s1=[]; decelMiss_s1=[]; decel_s1_phaseShift=[]; decel_s1_avgIBI=[];
%     accel_s0 = []; accelHit_s0=[]; accelMiss_s0=[]; accel_s0_phaseShift=[]; accel_s0_avgIBI=[];
%     decel_s0 = []; decelHit_s0=[]; decelMiss_s0=[]; decel_s0_phaseShift=[]; decel_s0_avgIBI=[];
%     accel_p1 = []; accelHit_p1=[]; accelMiss_p1=[]; accel_p1_phaseShift=[]; accel_p1_avgIBI=[];
%     decel_p1 = []; decelHit_p1=[]; decelMiss_p1=[]; decel_p1_phaseShift=[]; decel_p1_avgIBI=[];

    exhale_s2Phase = []; %these 6 values are for doing a phase locking analysis from only exhalation and only inhalation trials at the given time point!
    exhale_s1Phase = [];
    exhale_s0Phase = [];
    inhale_s2Phase = [];
    inhale_s1Phase = [];
    inhale_s0Phase = [];
    
    hit_brPeriods = [];
    miss_brPeriods = [];
    hitExhale_brPeriods = [];
    hitInhale_brPeriods = [];
    missExhale_brPeriods = [];
    missInhale_brPeriods = [];
    
    all_brPeriods = []; %this is the full breathing periods (inhale+exhale duration)
    hit_full_brs = [];
    miss_full_brs = [];
    
    %%% with the below variables, you are going to do circular-linear correlation plots for each subject!!!
    preStim_brPhases = []; preStim_brPhases_hit = []; preStim_brPhases_miss = []; preStim_brPhases_fast = []; preStim_brPhases_slow = [];
    preStim_ACDs = []; preStim_ACDs_hit = []; preStim_ACDs_miss = []; preStim_ACDs_fast = []; preStim_ACDs_slow = [];
    preStim_brPhases_highC = []; preStim_brPhases_lowC = []; preStim_ACDs_highC = []; preStim_ACDs_lowC = [];

    %%% with the below variables, you will separate the sine wave fits on s-2, s-1, and s-0 beats: to see the progression of the ACD!
    acd_s0_hit = []; acd_s0_miss = []; acd_s1_hit = []; acd_s1_miss = []; acd_s2_hit = []; acd_s2_miss = []; acd_s_p1_hit=[]; acd_s_p1_miss=[]; acd_s3_hit = []; acd_s3_miss = [];
    s_minZ_hitPhases=[]; s_minZ_missPhases=[]; s_min1_hitPhases=[]; s_min1_missPhases=[]; s_min2_hitPhases=[]; s_min2_missPhases=[]; s_min3_hitPhases=[]; s_min3_missPhases=[]; 
    s_p1_hitPhases=[]; s_p1_missPhases=[];

    %%% CONFIDENCE: with the below variables, you will separate the sine wave fits on s-2, s-1, and s-0 beats: to see the progression of the ACD!
    acd_s0_highC = []; acd_s0_lowC = []; acd_s1_highC = []; acd_s1_lowC = []; acd_s2_highC = []; acd_s2_lowC = []; acd_s_p1_highC=[]; acd_s_p1_lowC=[]; acd_s3_highC = []; acd_s3_lowC = [];
    s_minZ_highCPhases=[]; s_minZ_lowCPhases=[]; s_min1_highCPhases=[]; s_min1_lowCPhases=[]; s_min2_highCPhases=[]; s_min2_lowCPhases=[]; s_min3_highCPhases=[]; s_min3_lowCPhases=[]; 
    s_p1_highCPhases=[]; s_p1_lowCPhases=[];

    %%% 2x2 s0 vs s2 delta vertical offset analysis (HH: hit-highC, HL: hit-lowC, MH: miss-highC, ML: miss-lowC)
    acd_s0_HH = []; acd_s0_HL = []; acd_s0_MH = []; acd_s0_ML = [];
    acd_s2_HH = []; acd_s2_HL = []; acd_s2_MH = []; acd_s2_ML = [];
    s_minZ_HHPhases = []; s_minZ_HLPhases=[]; s_minZ_MHPhases=[]; s_minZ_MLPhases=[];
    s_min2_HHPhases = []; s_min2_HLPhases=[]; s_min2_MHPhases=[]; s_min2_MLPhases=[];
    
    fastResp_ibi_s2 = []; slowResp_ibi_s2 = [];
    fastResp_ibi_s1 = []; slowResp_ibi_s1 = [];
    fastResp_ibi_s0 = []; slowResp_ibi_s0 = [];
    fastResp_ibi_s_p1 = []; slowResp_ibi_s_p1 = [];

    highC_ibi_s2 = []; lowC_ibi_s2 = [];
    highC_ibi_s1 = []; lowC_ibi_s1 = [];
    highC_ibi_s0 = []; lowC_ibi_s0 = [];

    hitPLVs = [];
    missPLVs= [];

    hitRespTimes =[];
    missRespTimes=[];
    allRespTimes =[];

    goodLock_s2 = [];
    goodLockHit_s2 = [];
    goodLockMiss_s2 = [];
    badLock_s2 = [];
    badLockHit_s2 = [];
    badLockMiss_s2 = [];
    goodLock_s1 = [];
    goodLockHit_s1 = [];
    goodLockMiss_s1 = [];
    badLock_s1 = [];
    badLockHit_s1 = [];
    badLockMiss_s1 = [];
    goodLock_s0 = [];
    goodLockHit_s0 = [];
    goodLockMiss_s0 = [];
    badLock_s0 = [];
    badLockHit_s0 = [];
    badLockMiss_s0 = [];
    goodLock_s_p1 = [];
    goodLockHit_s_p1 = [];
    goodLockMiss_s_p1 = [];
    badLock_s_p1 = [];
    badLockHit_s_p1 = [];
    badLockMiss_s_p1 = [];
    
    
    addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\MATLAB\eeglab2023.0
    eeglab
    savePath = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\preliminary\VTD_EEG_preprocessing_newFilters\vtd' subjectNo]; %to save final HEP data
    directory = fullfile(savePath);

    fileList_br = dir(fullfile(directory, ['*breathing_contMerged_trueEvents_wTrough*.set']));

    breath_mergedCont = pop_loadset(fileList_br.name,savePath);
    breath_set = pop_epoch(breath_mergedCont,{'S  1'},[-6 10.001]);
    bad_events = contains({breath_set.event.type}, {'breathing peak'});
    breath_set = pop_editeventvals(breath_set,'delete', find(bad_events));
    %%%
    eventList = [];
    codeList = [];
    for i=1:size(breath_set.event,2)
        eventList{i} = breath_set.event(i).type; %list the event order (s1, rPeak, rPeak, s3, rPeak, s4...)
        codeList{i} = breath_set.event(i).code;
    end

    epochEvents = [];

    for i=1:size(breath_set.epoch,2) %scanning through each epoch
        currHit = nan;
        rPeakInds = [];
        preStimR = [];
        postStimR = [];
        for ev=1:size(breath_set.epoch(i).eventtype,2) %scanning through each event within the current epoch
            if strcmp(breath_set.epoch(i).eventtype{1,ev},'S  1') %because the event number should always be 1, as the epochs arte
                epochEvents{i,ev} = 'TS';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
%             elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'S  1') && ev~=1 %in subject 07, there are extra 'S  1' events, I am making them NULL. These are marker contaminations.
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'S  3') && (ev>=9 && ev<26)
                epochEvents{i,ev} = 'stim';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
                curr_stimOnset = breath_set.epoch(i).eventlatency{1,ev};
                stimInd = ev;
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'S  3') && (ev<9 || ev>=26)
                epochEvents{i,ev} = 'WrongStim';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'S  4')
                epochEvents{i,ev} = 'H';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'S  5')
                epochEvents{i,ev} = 'M';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'Sync Off')
                epochEvents{i,ev} = 'Sync Off';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'S  2') %there was one 'S  2' indicating the beginning of a stimulus-absent trial, but it is a marker contamination in this case.
                epochEvents{i,ev} = 'Null';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'r peak')
                rPeakInds = [rPeakInds ev];
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'resp_event')
                epochEvents{i,ev} = 'BRpeak';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            elseif strcmp(breath_set.epoch(i).eventtype{1,ev},'resp_trough')
                epochEvents{i,ev} = 'BRtrough';
                epochLatency{i,ev} = breath_set.epoch(i).eventlatency{1, ev};
            end
        end

        %% Here we order the r-peak events with respect to the stimulus onset event.

%         if i==310
%             stimInd =  7;
%             epochEvents{i,stimInd} = 'stim'; epochEvents{i,21} = 'WrongStim';
%         elseif i==208
%             stimInd = 8;
%             epochEvents{i,stimInd} = 'stim';
%         elseif i==335
%             epochEvents{i,9} = 'WrongStim';
%         end

%         if i==260
%             stimInd =  7;
%         end

%         if i==361
%             stimInd =  8;
%             epochEvents{i,stimInd} = 'stim'; epochEvents{i,10} = 'WrongStim'; epochEvents{i,23} = 'WrongStim';
%         end
        preStimR = rPeakInds(rPeakInds<stimInd);
        postStimR = rPeakInds(rPeakInds>stimInd);
        preCount=0;

        for c = numel(preStimR):-1:1
            if preCount==0
                epochEvents{i,preStimR(c)} = ['s' num2str(preCount)];
            else
                epochEvents{i,preStimR(c)} = ['s-' num2str(preCount)];
            end
            preCount = preCount+1;
        end

        postCount=0;
        for c=1:numel(postStimR)
            epochEvents{i,postStimR(c)} = ['s-p' num2str(c)];
        end

%         for ev=1:size(breath_set.epoch(i).eventtype,2) %scanning through each event within the current epoch
%             if strcmp(breath_set.epoch(i).eventtype{1,ev},'r peak')
%                 r_latency = breath_set.epoch(i).eventlatency{1,ev};
%                 if (r_latency <= curr_stimOnset-600) && (r_latency >= curr_stimOnset-1850)
%                     epochEvents{i,ev} = 'W';
%                 end
%             end
%         end

        for r=1:length(rPeakInds)-1 %to discard rPeak events that are followed by another rPeak in less than 600ms away!
            if breath_set.epoch(i).eventlatency{1,rPeakInds(r+1)}<breath_set.epoch(i).eventlatency{1,rPeakInds(r)}+600
                epochEvents{i,rPeakInds(r)} = 'Inv';
            end
        end

        %% testing hit vs miss for the current trial
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'stim')
                stim_ind = event;
                currHit=nan;
                for hm=(stim_ind+1):(stim_ind+7) %+7 is arbitrary atm
                    
                    if strcmp(epochEvents{i,hm},'H')
                        currHit = 1;
                        hitMiss(i) = 1;
                        
                        respLatency = epochLatency{i,hm};
                        respTime = respLatency - epochLatency{i,stim_ind};
                        hitRespTimes = [hitRespTimes respTime];
                        allRespTimes = [allRespTimes respTime];
                    elseif strcmp(epochEvents{i,hm},'M')
                        currHit = 0;
                        hitMiss(i) = 0;
                        respLatency = epochLatency{i,hm};
                        respTime = respLatency - epochLatency{i,stim_ind};
                        missRespTimes = [missRespTimes respTime];
                        allRespTimes = [allRespTimes respTime];
                    end
                end

                if i>1 && hitMiss(i-1)==1
                    prevHit=1;
                elseif i>1 && hitMiss(i-1)==0
                    prevHit=0;
                end

%                 if i==109
%                     currHit = 0;
%                     respLatency = epochLatency{i,34};
%                     respTime = respLatency - epochLatency{109,24};
%                     missRespTimes = [missRespTimes respTime];
%                     allRespTimes = [allRespTimes respTime];
%                 elseif i==348
%                     currHit = 1;
%                     respLatency = epochLatency{i,38};
%                     respTime = respLatency - epochLatency{348,29};
%                     hitRespTimes = [hitRespTimes respTime];
%                     allRespTimes = [allRespTimes respTime];
%                 end
                if isnan(currHit)
                    return
                end
            end
        end

        %% get breathing phase at the moment of stimulus onset: 07.02.25 --> getting this first because I may sort trials according to breathing phase at stimulus onset!
        br_ind=[];
        br_tr=[];
        br_trough=[];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'stim')
                targetBeat = event;
            elseif strcmp(epochEvents{i,event},'BRtrough')
                br_tr = [br_tr event];
            end
        end
        
        if ~isempty(targetBeat)
            br_pres_stim = br_ind(br_ind<targetBeat);
            br_pre = max(br_pres_stim);
    
            br_posts_stim = br_ind(br_ind>targetBeat);
            br_post = min(br_posts_stim);
            
            if ~isempty(br_pre) && ~isempty(br_post)
                br_trough = br_tr(br_tr>br_pre & br_tr<br_post);
            end

            if size(br_trough)>1
                br_trough=[];
            end
        else
            br_pre=[];
            br_post=[];
        end

        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_min_stimTime = epochLatency{i,targetBeat};
            br_period = post_brTime - pre_brTime;
            all_brPeriods = [all_brPeriods br_period];
    
            s_min_stimPhase = ((s_min_stimTime-pre_brTime)/(post_brTime-pre_brTime)) * 360;
    
            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_min_stimPhase= nan;
                br_exhalePeriod = nan;
                br_inhalePeriod = nan;
            elseif ~isempty(br_trough)
                brTroughTime = epochLatency{i,br_trough};
                br_exhalePeriod = brTroughTime - pre_brTime;
                br_inhalePeriod = post_brTime - brTroughTime;
                
                if currHit==1
                    hit_brPeriods = [hit_brPeriods br_exhalePeriod]; %br period is actually just exhale periods here
                    hit_full_brs = [hit_full_brs br_period]; %this is the full (inhale + exhale) breathing period !
                elseif currHit==0
                    miss_brPeriods = [miss_brPeriods br_exhalePeriod];
                    miss_full_brs = [miss_full_brs br_period];
                end
            end
            br_exhalePeriods(i,1) = br_exhalePeriod;
            br_inhalePeriods(i,1) = br_inhalePeriod;
            s_min_stim_allPhases(i,1) = s_min_stimPhase;
        else
            s_min_stimPhase=nan;
            s_min_stim_allPhases(i,1) =nan;
        end
        
        if s_min_stimPhase>0 && s_min_stimPhase<=180
            if currHit == 1
                exhaleHit_stim = [exhaleHit_stim i];
            elseif currHit == 0
                exhaleMiss_stim = [exhaleMiss_stim i];
            end
            if currHit ==1 && ~isempty(br_trough) %to compare the average breathing periods in hit vs miss trials.
                hitExhale_brPeriods = [hitExhale_brPeriods br_exhalePeriod];
            elseif currHit ==0
                missExhale_brPeriods = [missExhale_brPeriods br_exhalePeriod];
            end
        elseif s_min_stimPhase>180 && s_min_stimPhase<=360
            if currHit == 1
                inhaleHit_stim = [inhaleHit_stim i];
            elseif currHit == 0
                inhaleMiss_stim = [inhaleMiss_stim i];
            end
            if currHit ==1 %to compare the average breathing periods in hit vs miss trials.
                hitInhale_brPeriods = [hitInhale_brPeriods br_exhalePeriod];
            elseif currHit ==0
                missInhale_brPeriods = [missInhale_brPeriods br_exhalePeriod];
            end
        end

        %% get breathing phase at the moment of S-3 r-peak
        br_ind=[];
        previousBeat=[];
        targetBeat = [];
        nextBeat = [];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'s-3')
                targetBeat = event;
            elseif strcmp(epochEvents{i,event},'s-4')
                previousBeat = event;
            elseif strcmp(epochEvents{i,event},'s-2')
                nextBeat = event;
            end
        end
        if ~isempty(previousBeat) && ~isempty(nextBeat) && ~isempty(targetBeat)
            ibi_prev = epochLatency{i,targetBeat} - epochLatency{i,previousBeat};
            ibi_s3(i) = ibi_prev;
            ibi_next = epochLatency{i,nextBeat} - epochLatency{i,targetBeat};
            diff_s3 = ibi_next - ibi_prev;
            acd_s3(i) = diff_s3;
        else
            ibi_s3(i) = nan;
            acd_s3(i)=nan;
            diff_s3=nan;
        end
        
        if ~isempty(targetBeat)
            br_pres3 = br_ind(br_ind<targetBeat);
            br_pre = max(br_pres3);
    
            br_posts3 = br_ind(br_ind>targetBeat);
            br_post = min(br_posts3);
        else
            br_pre = [];
            br_post = [];
        end
        
        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_min3Time = epochLatency{i,targetBeat};

            s_min3Phase = ((s_min3Time-pre_brTime)/(post_brTime-pre_brTime)) * 360;

            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_min3Phase=nan;
            end

            s_min3_allPhases(i,1) = s_min3Phase;
        else
            s_min3Phase=nan;
            s_min3_allPhases(i,1) = nan;
        end
        preStim_brPhases = [preStim_brPhases s_min3Phase];
        preStim_ACDs = [preStim_ACDs diff_s3];
        if currHit==1
            preStim_brPhases_hit = [preStim_brPhases_hit s_min3Phase];
            preStim_ACDs_hit = [preStim_ACDs_hit diff_s3];
            acd_s3_hit = [acd_s3_hit diff_s3];
            s_min3_hitPhases = [s_min3_hitPhases s_min3Phase];
        elseif currHit==0
            preStim_brPhases_miss = [preStim_brPhases_miss s_min3Phase];
            preStim_ACDs_miss = [preStim_ACDs_miss diff_s3];
            acd_s3_miss = [acd_s3_miss diff_s3];
            s_min3_missPhases = [s_min3_missPhases s_min3Phase];
        end

        %% get breathing phase at the moment of S-2 r-peak
        br_ind=[];
        previousBeat=[];
        targetBeat=[];
        nextBeat=[];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'s-2')
                targetBeat = event;
            elseif strcmp(epochEvents{i,event},'s-3')
                previousBeat = event;
            elseif strcmp(epochEvents{i,event},'s-1')
                nextBeat = event;
            end
        end
        
        if ~isempty(previousBeat) && ~isempty(nextBeat) && ~isempty(targetBeat)
            ibi_prev = epochLatency{i,targetBeat} - epochLatency{i,previousBeat};
            ibi_s2(i) = ibi_prev;
            ibi_next = epochLatency{i,nextBeat} - epochLatency{i,targetBeat};
            diff_s2 = ibi_next - ibi_prev;
            acd_s2(i) = diff_s2;
        else
            ibi_s2(i) = nan;
            acd_s2(i)=nan;
            diff_s2=nan;
        end
        
        if ~isempty(targetBeat)
            br_pres2 = br_ind(br_ind<targetBeat);
            br_pre = max(br_pres2);
    
            br_posts2 = br_ind(br_ind>targetBeat);
            br_post = min(br_posts2);
        else
            br_pre = [];
            br_post = [];
        end
        
        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_min2Time = epochLatency{i,targetBeat};

            s_min2Phase = ((s_min2Time-pre_brTime)/(post_brTime-pre_brTime)) * 360;

            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_min2Phase=nan;
            end

            s_min2_allPhases(i,1) = s_min2Phase;
        else
            s_min2Phase=nan;
            s_min2_allPhases(i,1) =nan;
        end
        preStim_brPhases = [preStim_brPhases s_min2Phase];
        preStim_ACDs = [preStim_ACDs diff_s2];
        if currHit==1
            preStim_brPhases_hit = [preStim_brPhases_hit s_min2Phase];
            preStim_ACDs_hit = [preStim_ACDs_hit diff_s2];
            acd_s2_hit = [acd_s2_hit diff_s2];
            s_min2_hitPhases = [s_min2_hitPhases s_min2Phase];
            if ismember(i,highC_trs)
                acd_s2_HH = [acd_s2_HH diff_s2];
                s_min2_HHPhases = [s_min2_HHPhases s_min2Phase];
            elseif ismember(i,lowC_trs)
                acd_s2_HL = [acd_s2_HL diff_s2];
                s_min2_HLPhases = [s_min2_HLPhases s_min2Phase];
            end
        elseif currHit==0
            preStim_brPhases_miss = [preStim_brPhases_miss s_min2Phase];
            preStim_ACDs_miss = [preStim_ACDs_miss diff_s2];
            acd_s2_miss = [acd_s2_miss diff_s2];
            s_min2_missPhases = [s_min2_missPhases s_min2Phase];
            if ismember(i,highC_trs)
                acd_s2_MH = [acd_s2_MH diff_s2];
                s_min2_MHPhases = [s_min2_MHPhases s_min2Phase];
            elseif ismember(i,lowC_trs)
                acd_s2_ML = [acd_s2_ML diff_s2];
                s_min2_MLPhases = [s_min2_MLPhases s_min2Phase];
            end
        end
        if ismember(i,resps.fastResps)
            fastResp_ibi_s2 = [fastResp_ibi_s2 ibi_s2(i)];
            preStim_brPhases_fast = [preStim_brPhases_fast s_min2Phase];
            preStim_ACDs_fast = [preStim_ACDs_fast diff_s2];
        elseif ismember(i,resps.slowResps)
            slowResp_ibi_s2 = [slowResp_ibi_s2 ibi_s2(i)];
            preStim_brPhases_slow = [preStim_brPhases_slow s_min2Phase];
            preStim_ACDs_slow = [preStim_ACDs_slow diff_s2];
        end

        if ismember(i,highC_trs)
            highC_ibi_s2 = [highC_ibi_s2 ibi_s2(i)];
            preStim_brPhases_highC = [preStim_brPhases_highC s_min2Phase];
            preStim_ACDs_highC = [preStim_ACDs_highC diff_s2];
            acd_s2_highC = [acd_s2_highC diff_s2];
            s_min2_highCPhases = [s_min2_highCPhases s_min2Phase];
        elseif ismember(i,lowC_trs)
            lowC_ibi_s2 = [lowC_ibi_s2 ibi_s2(i)];
            preStim_brPhases_lowC = [preStim_brPhases_lowC s_min2Phase];
            preStim_ACDs_lowC = [preStim_ACDs_lowC diff_s2];
            acd_s2_lowC = [acd_s2_lowC diff_s2];
            s_min2_lowCPhases = [s_min2_lowCPhases s_min2Phase];
        end


        if s_min2Phase>0 && s_min2Phase<=180
            exhale_s2 = [exhale_s2 i];
            exhale_s2Phase = [exhale_s2Phase s_min2Phase];
            if currHit ==1
                exhaleHit_s2 = [exhaleHit_s2 i];
            elseif currHit ==0
                exhaleMiss_s2 = [exhaleMiss_s2 i];
            end
        elseif s_min2Phase>180 && s_min2Phase<=360
            inhale_s2 = [inhale_s2 i];
            inhale_s2Phase = [inhale_s2Phase s_min2Phase];
            if currHit ==1
                inhaleHit_s2 = [inhaleHit_s2 i];
            elseif currHit ==0
                inhaleMiss_s2 = [inhaleMiss_s2 i];
            end
        end

%         if (deg2rad(s_min2Phase)>7*pi/6 && deg2rad(s_min3Phase)>7*pi/6)
%             if deg2rad(s_min2Phase)>deg2rad(s_min3Phase) %acceleration is expected in this case (advancing within inhalation)
%                 accel_s2 = [accel_s2 i];
%                 accel_s2_phaseShift = [accel_s2_phaseShift (s_min2Phase - s_min3Phase)]; %phase shift will tell us if there is more phase jump in expected deceleration than acceleration; which could explain the extra cardiac deceleration
%                 accel_s2_avgIBI = [accel_s2_avgIBI ibi_prev];
%                 if currHit ==1
%                     accelHit_s2 = [accelHit_s2 i];
%                 elseif currHit ==0
%                     accelMiss_s2 = [accelMiss_s2 i];
%                 end
%             end
%         elseif (deg2rad(s_min2Phase)>0 && deg2rad(s_min3Phase)>0) && (deg2rad(s_min2Phase)<5*pi/6 && deg2rad(s_min3Phase)<5*pi/6) 
%             if deg2rad(s_min2Phase)>deg2rad(s_min3Phase) %deceleration is expected in this case (advancing within exhalation)
%                 decel_s2 = [decel_s2 i];
%                 decel_s2_phaseShift = [decel_s2_phaseShift (s_min2Phase - s_min3Phase)];
%                 decel_s2_avgIBI = [decel_s2_avgIBI ibi_prev];
%                 if currHit ==1
%                     decelHit_s2 = [decelHit_s2 i];
%                 elseif currHit ==0
%                     decelMiss_s2 = [decelMiss_s2 i];
%                 end
%             end
%         end
% 
%         if s_min_stimPhase>90 && s_min_stimPhase<=225 %a 135 degrees of phase, which is definitely a good phase locking
%             goodLock_s2 = [goodLock_s2 i];
%             if currHit ==1
%                 goodLockHit_s2 = [goodLockHit_s2 i];
%             elseif currHit ==0
%                 goodLockMiss_s2 = [goodLockMiss_s2 i];
%             end
%         elseif s_min_stimPhase<45 || s_min_stimPhase>270 %a 135 degrees of phase, which is definitely a bad phase locking
%             badLock_s2 = [badLock_s2 i];
%             if currHit ==1
%                 badLockHit_s2 = [badLockHit_s2 i];
%             elseif currHit ==0
%                 badLockMiss_s2 = [badLockMiss_s2 i];
%             end
%         end

        %% get breathing phase at the moment of S-1 r-peak
        br_ind=[];
        previousBeat=[];
        targetBeat=[];
        nextBeat=[];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'s-1')
                targetBeat = event;
            elseif strcmp(epochEvents{i,event},'s0')
                nextBeat = event;
            elseif strcmp(epochEvents{i,event},'s-2')
                previousBeat = event;
            end
        end
        
        if ~isempty(previousBeat) && ~isempty(nextBeat) && ~isempty(targetBeat)
            ibi_prev = epochLatency{i,targetBeat} - epochLatency{i,previousBeat};
            ibi_next = epochLatency{i,nextBeat} - epochLatency{i,targetBeat};
            ibi_s1(i) = ibi_prev;
            diff_s1 = ibi_next - ibi_prev;
            acd_s1(i) = diff_s1;
        else
            ibi_s1(i) = nan;
            acd_s1(i) = nan;
            diff_s1 = nan;
        end
        
        if ~isempty(targetBeat)
            br_pres1 = br_ind(br_ind<targetBeat);
            br_pre = max(br_pres1);
    
            br_posts1 = br_ind(br_ind>targetBeat);
            br_post = min(br_posts1);
        else
            br_pre = [];
            br_post = [];
        end

        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_min1Time = epochLatency{i,targetBeat};
    
            s_min1Phase = ((s_min1Time-pre_brTime)/(post_brTime-pre_brTime)) * 360;
    
            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_min1Phase=nan;
            end

            s_min1_allPhases(i,1) = s_min1Phase;
        else
            s_min1Phase=nan;
            s_min1_allPhases(i,1) =nan;
        end
        preStim_brPhases = [preStim_brPhases s_min1Phase];
        preStim_ACDs = [preStim_ACDs diff_s1];
        if currHit==1
            preStim_brPhases_hit = [preStim_brPhases_hit s_min1Phase];
            preStim_ACDs_hit = [preStim_ACDs_hit diff_s1];
            acd_s1_hit = [acd_s1_hit diff_s1];
            s_min1_hitPhases = [s_min1_hitPhases s_min1Phase];
        elseif currHit==0
            preStim_brPhases_miss = [preStim_brPhases_miss s_min1Phase];
            preStim_ACDs_miss = [preStim_ACDs_miss diff_s1];
            acd_s1_miss = [acd_s1_miss diff_s1];
            s_min1_missPhases = [s_min1_missPhases s_min1Phase];
        end
        if ismember(i,resps.fastResps)
            fastResp_ibi_s1 = [fastResp_ibi_s1 ibi_s1(i)];
            preStim_brPhases_fast = [preStim_brPhases_fast s_min1Phase];
            preStim_ACDs_fast = [preStim_ACDs_fast diff_s1];
        elseif ismember(i,resps.slowResps)
            slowResp_ibi_s1 = [slowResp_ibi_s1 ibi_s1(i)];
            preStim_brPhases_slow = [preStim_brPhases_slow s_min1Phase];
            preStim_ACDs_slow = [preStim_ACDs_slow diff_s1];
        end
        if ismember(i,highC_trs)
            highC_ibi_s1 = [highC_ibi_s1 ibi_s1(i)];
            preStim_brPhases_highC = [preStim_brPhases_highC s_min1Phase];
            preStim_ACDs_highC = [preStim_ACDs_highC diff_s1];
            acd_s1_highC = [acd_s1_highC diff_s1];
            s_min1_highCPhases = [s_min1_highCPhases s_min1Phase];
        elseif ismember(i,lowC_trs)
            lowC_ibi_s1 = [lowC_ibi_s1 ibi_s1(i)];
            preStim_brPhases_lowC = [preStim_brPhases_lowC s_min1Phase];
            preStim_ACDs_lowC = [preStim_ACDs_lowC diff_s1];
            acd_s1_lowC = [acd_s1_lowC diff_s1];
            s_min1_lowCPhases = [s_min1_lowCPhases s_min1Phase];
        end

        if s_min1Phase>0 && s_min1Phase<=180
            exhale_s1 = [exhale_s1 i];
            exhale_s1Phase = [exhale_s1Phase s_min1Phase];
            if currHit ==1
                exhaleHit_s1 = [exhaleHit_s1 i];
            elseif currHit ==0
                exhaleMiss_s1 = [exhaleMiss_s1 i];
            end
        elseif s_min1Phase>180 && s_min1Phase<=360
            inhale_s1 = [inhale_s1 i];
            inhale_s1Phase = [inhale_s1Phase s_min1Phase];
            if currHit ==1
                inhaleHit_s1 = [inhaleHit_s1 i];
            elseif currHit ==0
                inhaleMiss_s1 = [inhaleMiss_s1 i];
            end
        end

        %%%% 16.06.25 addition: if rad(angle) of min1 is bigger than that of min2, this means that cardiac deceleration is favored in this trial via
        %%%% RSA
        %Nested ifs take care of the following: if there is progress either within inhalation, or within exhalation...
%         if (deg2rad(s_min1Phase)>7*pi/6 && deg2rad(s_min2Phase)>7*pi/6)
%             if deg2rad(s_min1Phase)>deg2rad(s_min2Phase) %acceleration is expected in this case (advancing within inhalation)
%                 accel_s1 = [accel_s1 i];
%                 accel_s1_phaseShift = [accel_s1_phaseShift (s_min1Phase-s_min2Phase)];
%                 accel_s1_avgIBI = [accel_s1_avgIBI ibi_prev];
%                 if currHit ==1
%                     accelHit_s1 = [accelHit_s1 i];
%                 elseif currHit ==0
%                     accelMiss_s1 = [accelMiss_s1 i];
%                 end
%             end
%         elseif (deg2rad(s_min1Phase)>0 && deg2rad(s_min2Phase)>0) && (deg2rad(s_min1Phase)<5*pi/6 && deg2rad(s_min2Phase)<5*pi/6)
%             if deg2rad(s_min1Phase)>deg2rad(s_min2Phase) %deceleration is expected in this case (advancing within exhalation)
%                 decel_s1 = [decel_s1 i];
%                 decel_s1_phaseShift = [decel_s1_phaseShift (s_min1Phase-s_min2Phase)];
%                 decel_s1_avgIBI = [decel_s1_avgIBI ibi_prev];
%                 if currHit ==1
%                     decelHit_s1 = [decelHit_s1 i];
%                 elseif currHit ==0
%                     decelMiss_s1 = [decelMiss_s1 i];
%                 end
%             end
%         end
% 
%         if s_min_stimPhase>90 && s_min_stimPhase<=225 %a 135 degrees of phase, which is definitely a good phase locking
%             goodLock_s1 = [goodLock_s1 i];
%             if currHit ==1
%                 goodLockHit_s1 = [goodLockHit_s1 i];
%             elseif currHit ==0
%                 goodLockMiss_s1 = [goodLockMiss_s1 i];
%             end
%         elseif s_min_stimPhase<45 || s_min_stimPhase>270 %a 135 degrees of phase, which is definitely a bad phase locking
%             badLock_s1 = [badLock_s1 i];
%             if currHit ==1
%                 badLockHit_s1 = [badLockHit_s1 i];
%             elseif currHit ==0
%                 badLockMiss_s1 = [badLockMiss_s1 i];
%             end
%         end

        %% get breathing phase at the moment of S-0
        br_ind=[];
        previousBeat=[];
        targetBeat=[];
        nextBeat=[];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'s0')
                targetBeat = event;
            elseif strcmp(epochEvents{i,event},'s-p1')
                nextBeat = event;
            elseif strcmp(epochEvents{i,event},'s-1')
                previousBeat = event;
            end
        end
        
        if ~isempty(previousBeat) && ~isempty(nextBeat) && ~isempty(targetBeat)
            ibi_prev = epochLatency{i,targetBeat} - epochLatency{i,previousBeat};
            ibi_next = epochLatency{i,nextBeat} - epochLatency{i,targetBeat};
            ibi_s0(i) = ibi_prev;
            diff_s0 = ibi_next - ibi_prev;
            acd_s0(i) = diff_s0;
        else
            ibi_s0(i) = nan;
            acd_s0(i) = nan;
            diff_s0 = nan;
        end
        
        if ~isempty(targetBeat)
            br_presZ = br_ind(br_ind<targetBeat);
            br_pre = max(br_presZ);
    
            br_postsZ = br_ind(br_ind>targetBeat);
            br_post = min(br_postsZ);
        else
            br_pre =[];
            br_post=[];
        end

        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_minZTime = epochLatency{i,targetBeat};

            s_minZPhase = ((s_minZTime-pre_brTime)/(post_brTime-pre_brTime)) * 360;

            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_minZPhase=nan;
            end

            s_minZ_allPhases(i,1) = s_minZPhase;
        else
            s_minZPhase=nan;
            s_minZ_allPhases(i,1) =nan;
        end
        preStim_brPhases = [preStim_brPhases s_minZPhase];
        preStim_ACDs = [preStim_ACDs diff_s0];
        if currHit==1
            preStim_brPhases_hit = [preStim_brPhases_hit s_minZPhase];
            preStim_ACDs_hit = [preStim_ACDs_hit diff_s0];
            acd_s0_hit = [acd_s0_hit diff_s0];
            s_minZ_hitPhases = [s_minZ_hitPhases s_minZPhase];
            if ismember(i,highC_trs)
                acd_s0_HH = [acd_s0_HH diff_s0];
                s_minZ_HHPhases = [s_minZ_HHPhases s_minZPhase];
            elseif ismember(i,lowC_trs)
                acd_s0_HL = [acd_s0_HL diff_s0];
                s_minZ_HLPhases = [s_minZ_HLPhases s_minZPhase];
            end
        elseif currHit==0
            preStim_brPhases_miss = [preStim_brPhases_miss s_minZPhase];
            preStim_ACDs_miss = [preStim_ACDs_miss diff_s0];
            acd_s0_miss = [acd_s0_miss diff_s0];
            s_minZ_missPhases = [s_minZ_missPhases s_minZPhase];
            if ismember(i,highC_trs)
                acd_s0_MH = [acd_s0_MH diff_s0];
                s_minZ_MHPhases = [s_minZ_MHPhases s_minZPhase];
            elseif ismember(i,lowC_trs)
                acd_s0_ML = [acd_s0_ML diff_s0];
                s_minZ_MLPhases = [s_minZ_MLPhases s_minZPhase];
            end
        end
        if ismember(i,resps.fastResps)
            fastResp_ibi_s0 = [fastResp_ibi_s0 ibi_s0(i)];
            preStim_brPhases_fast = [preStim_brPhases_fast s_minZPhase];
            preStim_ACDs_fast = [preStim_ACDs_fast diff_s0];
        elseif ismember(i,resps.slowResps)
            slowResp_ibi_s0 = [slowResp_ibi_s0 ibi_s0(i)];
            preStim_brPhases_slow = [preStim_brPhases_slow s_minZPhase];
            preStim_ACDs_slow = [preStim_ACDs_slow diff_s0];
        end
        if ismember(i,highC_trs)
            highC_ibi_s0 = [highC_ibi_s0 ibi_s0(i)];
            preStim_brPhases_highC = [preStim_brPhases_highC s_minZPhase];
            preStim_ACDs_highC = [preStim_ACDs_highC diff_s0];
            acd_s0_highC = [acd_s0_highC diff_s0];
            s_minZ_highCPhases = [s_minZ_highCPhases s_minZPhase];
        elseif ismember(i,lowC_trs)
            lowC_ibi_s0 = [lowC_ibi_s0 ibi_s0(i)];
            preStim_brPhases_lowC = [preStim_brPhases_lowC s_minZPhase];
            preStim_ACDs_lowC = [preStim_ACDs_lowC diff_s0 ];
            acd_s0_lowC = [acd_s0_lowC diff_s0];
            s_minZ_lowCPhases = [s_minZ_lowCPhases s_minZPhase];
        end


        if s_minZPhase>0 && s_minZPhase<=180
            exhale_s0 = [exhale_s0 i];
            exhale_s0Phase = [exhale_s0Phase s_minZPhase];
            if currHit ==1
                exhaleHit_s0 = [exhaleHit_s0 i];
            elseif currHit ==0
                exhaleMiss_s0 = [exhaleMiss_s0 i];
            end
        elseif s_minZPhase>180 && s_minZPhase<=360
            inhale_s0 = [inhale_s0 i];
            inhale_s0Phase = [inhale_s0Phase s_minZPhase];
            if currHit ==1
                inhaleHit_s0 = [inhaleHit_s0 i];
            elseif currHit ==0
                inhaleMiss_s0 = [inhaleMiss_s0 i];
            end
        end

%         if (deg2rad(s_minZPhase)>7*pi/6 && deg2rad(s_min1Phase)>7*pi/6)
%             if deg2rad(s_minZPhase)>deg2rad(s_min1Phase) %acceleration is expected in this case (advancing within inhalation)
%                 accel_s0 = [accel_s0 i];
%                 accel_s0_phaseShift = [accel_s0_phaseShift (s_minZPhase-s_min1Phase)];
%                 accel_s0_avgIBI = [accel_s0_avgIBI ibi_prev];
%                 if currHit ==1
%                     accelHit_s0 = [accelHit_s0 i];
%                 elseif currHit ==0
%                     accelMiss_s0 = [accelMiss_s0 i];
%                 end
%             end
%         elseif (deg2rad(s_minZPhase)>0 && deg2rad(s_min1Phase)>0) && (deg2rad(s_minZPhase)<5*pi/6 && deg2rad(s_min1Phase)<5*pi/6)
%             if deg2rad(s_minZPhase)>deg2rad(s_min1Phase) %deceleration is expected in this case (advancing within exhalation)
%                 decel_s0 = [decel_s0 i];
%                 decel_s0_phaseShift = [decel_s0_phaseShift (s_minZPhase-s_min1Phase)];
%                 decel_s0_avgIBI = [decel_s0_avgIBI ibi_prev];
%                 if currHit ==1
%                     decelHit_s0 = [decelHit_s0 i];
%                 elseif currHit ==0
%                     decelMiss_s0 = [decelMiss_s0 i];
%                 end
%             end
%         end
% 
%         if s_min_stimPhase>90 && s_min_stimPhase<=225 %a 135 degrees of phase, which is definitely a good phase locking
%             goodLock_s0 = [goodLock_s0 i];
%             if currHit ==1
%                 goodLockHit_s0 = [goodLockHit_s0 i];
%             elseif currHit ==0
%                 goodLockMiss_s0 = [goodLockMiss_s0 i];
%             end
%         elseif s_min_stimPhase<45 || s_min_stimPhase>270 %a 135 degrees of phase, which is definitely a bad phase locking
%             badLock_s0 = [badLock_s0 i];
%             if currHit ==1
%                 badLockHit_s0 = [badLockHit_s0 i];
%             elseif currHit ==0
%                 badLockMiss_s0 = [badLockMiss_s0 i];
%             end
%         end


        %% get breathing phase at the moment of S-plus 1 r-peak
        br_ind=[];
        previousBeat=[];
        targetBeat=[];
        nextBeat=[];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'s-p1')
                targetBeat = event;
            elseif strcmp(epochEvents{i,event},'s-p2')
                nextBeat = event;
            elseif strcmp(epochEvents{i,event},'s0')
                previousBeat = event;
            end
        end
        
        if ~isempty(previousBeat) && ~isempty(nextBeat) && ~isempty(targetBeat)
            ibi_prev = epochLatency{i,targetBeat} - epochLatency{i,previousBeat};
            ibi_next = epochLatency{i,nextBeat} - epochLatency{i,targetBeat};
            ibi_s_p1(i) = ibi_prev;
            diff_s_p1 = ibi_next - ibi_prev;
            acd_s_p1(i) = diff_s_p1;
        else
            ibi_s_p1(i) = nan;
            acd_s_p1(i) = nan;
        end
        
        if ~isempty(targetBeat)
            br_pres_p1 = br_ind(br_ind<targetBeat);
            br_pre = max(br_pres_p1);
    
            br_posts_p1 = br_ind(br_ind>targetBeat);
            br_post = min(br_posts_p1);
        else
            br_pre = [];
            br_post = [];
        end

        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_plus1Time = epochLatency{i,targetBeat};
    
            s_plus1Phase = ((s_plus1Time-pre_brTime)/(post_brTime-pre_brTime)) * 360;
    
            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_plus1Phase=nan;
            end

            s_plus1_allPhases(i,1) = s_plus1Phase;
        else
            s_plus1Phase=nan;
            s_plus1_allPhases(i,1) =nan;
        end

        if currHit==1
            acd_s_p1_hit = [acd_s_p1_hit diff_s_p1];
            s_p1_hitPhases = [s_p1_hitPhases s_plus1Phase];
        elseif currHit==0
            acd_s_p1_miss = [acd_s_p1_miss diff_s_p1];
            s_p1_missPhases = [s_p1_missPhases s_plus1Phase];
        end

        if ismember(i,resps.fastResps)
            fastResp_ibi_s_p1 = [fastResp_ibi_s_p1 ibi_s_p1(i)];
        elseif ismember(i,resps.slowResps)
            slowResp_ibi_s_p1 = [slowResp_ibi_s_p1 ibi_s_p1(i)];
        end

        if s_plus1Phase>0 && s_plus1Phase<=180
            exhale_s_p1 = [exhale_s_p1 i];
            if currHit ==1
                exhaleHit_s_p1 = [exhaleHit_s_p1 i];
            elseif currHit ==0
                exhaleMiss_s_p1 = [exhaleMiss_s_p1 i];
            end
        elseif s_plus1Phase>180 && s_plus1Phase<=360
            inhale_s_p1 = [inhale_s_p1 i];
            if currHit ==1
                inhaleHit_s_p1 = [inhaleHit_s_p1 i];
            elseif currHit ==0
                inhaleMiss_s_p1 = [inhaleMiss_s_p1 i];
            end
        end

%         if (deg2rad(s_plus1Phase)>7*pi/6 && deg2rad(s_minZPhase)>7*pi/6)
%             if deg2rad(s_plus1Phase)>deg2rad(s_minZPhase) %deceleration is expected in this case (advancing within inhalation)
%                 accel_p1 = [accel_p1 i];
%                 accel_p1_phaseShift = [accel_p1_phaseShift (s_plus1Phase-s_minZPhase)];
%                 accel_p1_avgIBI = [accel_p1_avgIBI ibi_prev];
%                 if currHit ==1
%                     accelHit_p1 = [accelHit_p1 i];
%                 elseif currHit ==0
%                     accelMiss_p1 = [accelMiss_p1 i];
%                 end
%             end
%         elseif (deg2rad(s_plus1Phase)>0 && deg2rad(s_minZPhase)>0) && (deg2rad(s_plus1Phase)<5*pi/6 && deg2rad(s_minZPhase)<5*pi/6)
%             if deg2rad(s_plus1Phase)>deg2rad(s_minZPhase) %deceleration is expected in this case (advancing within exhalation)
%                 decel_p1 = [decel_p1 i];
%                 decel_p1_phaseShift = [decel_p1_phaseShift (s_plus1Phase-s_minZPhase)];
%                 decel_p1_avgIBI = [decel_p1_avgIBI ibi_prev];
%                 if currHit ==1
%                     decelHit_p1 = [decelHit_p1 i];
%                 elseif currHit ==0
%                     decelMiss_p1 = [decelMiss_p1 i];
%                 end
%             end
%         end
% 
%         if s_min_stimPhase>90 && s_min_stimPhase<=225 %a 135 degrees of phase, which is definitely a good phase locking
%             goodLock_s_p1 = [goodLock_s_p1 i];
%             if currHit ==1
%                 goodLockHit_s_p1 = [goodLockHit_s_p1 i];
%             elseif currHit ==0
%                 goodLockMiss_s_p1 = [goodLockMiss_s_p1 i];
%             end
%         elseif s_min_stimPhase<45 || s_min_stimPhase>270 %a 135 degrees of phase, which is definitely a bad phase locking
%             badLock_s_p1 = [badLock_s_p1 i];
%             if currHit ==1
%                 badLockHit_s_p1 = [badLockHit_s_p1 i];
%             elseif currHit ==0
%                 badLockMiss_s_p1 = [badLockMiss_s_p1 i];
%             end
%         end

        %% get breathing phase at the moment of S-plus 2 r-peak
        br_ind=[];
        for event=1:size(epochEvents,2)
            if strcmp(epochEvents{i,event},'BRpeak')
                br_ind = [br_ind event];
            elseif strcmp(epochEvents{i,event},'s-p2')
                targetBeat = event;
            end
        end
        
        if ~isempty(targetBeat)
            br_pres_p2 = br_ind(br_ind<targetBeat);
            br_pre = max(br_pres_p2);
    
            br_posts_p2 = br_ind(br_ind>targetBeat);
            br_post = min(br_posts_p2);
        else
            br_pre = [];
            br_post = [];
        end

        if ~isempty(br_pre) && ~isempty(br_post)
            pre_brTime = epochLatency{i,br_pre};
            post_brTime = epochLatency{i,br_post};

            s_plus2Time = epochLatency{i,targetBeat};
    
            s_plus2Phase = ((s_plus2Time-pre_brTime)/(post_brTime-pre_brTime)) * 360;
    
            if ((post_brTime - pre_brTime) > 7000) || ((post_brTime - pre_brTime) < 1500)
                s_plus2Phase=nan;
            end

            s_plus2_allPhases(i,1) = s_plus2Phase;
        else
            s_plus2Phase=nan;
            s_plus2_allPhases(i,1) =nan;
        end

        clearvars peaks r_peaks PLV times phases resp_phase_at_rpeaks phase_differences midpoints midpoint_phases currEvents currLatency

    end
    
    eventListNew = reshape(epochEvents',1,size(epochEvents,1)*size(epochEvents,2));
    eventListCoded = eventListNew(~cellfun(@isempty, eventListNew));

    %%% Check that eventListCoded has the same size as eventList, making sure that each event is coded according to your new scheme.
    if size(eventList,2) ~= size(eventListCoded,2)
        disp(['Your new event codes are erroneous: Not all events are coded'])
        return
    end
    
    %%% Load this new version of event coding into your EEGLAB dataset, to visualize and play around with ERPs easier!
    for i=1:size(breath_set.event,2)
        breath_set.event(i).type = eventListCoded{i}; %put the new event order into the EEG.event struct, to be read by EEGLAB
    end
    
    %% This part gets the ACD values at the corresponding time period, in the trials where the subject is exhaling or inhaling...

    acd_exhale_s2 = acd_s2(exhale_s2);
    acd_exhaleHit_s2 = acd_s2(exhaleHit_s2);
    acd_exhaleMiss_s2 = acd_s2(exhaleMiss_s2);
    acd_exhale_s1 = acd_s1(exhale_s1);
    acd_exhaleHit_s1 = acd_s1(exhaleHit_s1);
    acd_exhaleMiss_s1 = acd_s1(exhaleMiss_s1);
    acd_exhale_s0 = acd_s0(exhale_s0);
    acd_exhaleHit_s0 = acd_s0(exhaleHit_s0);
    acd_exhaleMiss_s0 = acd_s0(exhaleMiss_s0);
    acd_exhale_s_p1 = acd_s_p1(exhale_s_p1);
    acd_exhaleHit_s_p1 = acd_s_p1(exhaleHit_s_p1);
    acd_exhaleMiss_s_p1 = acd_s_p1(exhaleMiss_s_p1);

    exh_s2_hitRate = 100*size(exhaleHit_s2,2)/(size(exhaleHit_s2,2) + size(exhaleMiss_s2,2));
    inh_s2_hitRate = 100*size(inhaleHit_s2,2)/(size(inhaleHit_s2,2) + size(inhaleMiss_s2,2));
    exh_s1_hitRate = 100*size(exhaleHit_s1,2)/(size(exhaleHit_s1,2) + size(exhaleMiss_s1,2));
    inh_s1_hitRate = 100*size(inhaleHit_s1,2)/(size(inhaleHit_s1,2) + size(inhaleMiss_s1,2));
    exh_s0_hitRate = 100*size(exhaleHit_s0,2)/(size(exhaleHit_s0,2) + size(exhaleMiss_s0,2));
    inh_s0_hitRate = 100*size(inhaleHit_s0,2)/(size(inhaleHit_s0,2) + size(inhaleMiss_s0,2));
    exh_stim_hitRate = 100*size(exhaleHit_stim,2)/(size(exhaleHit_stim,2) + size(exhaleMiss_stim,2));
    inh_stim_hitRate = 100*size(inhaleHit_stim,2)/(size(inhaleHit_stim,2) + size(inhaleMiss_stim,2));
    
    mean_hitBrPeriods = mean(hit_brPeriods,'omitnan');
    mean_missBrPeriods = mean(miss_brPeriods,'omitnan');
    mean_hitExhaleBrPeriods = mean(hitExhale_brPeriods,'omitnan');
    mean_hitInhaleBrPeriods = mean(hitInhale_brPeriods,'omitnan');
    mean_missExhaleBrPeriods = mean(missExhale_brPeriods,'omitnan');
    mean_missInhaleBrPeriods = mean(missInhale_brPeriods,'omitnan');

    mean_all_brs = mean(all_brPeriods,'omitnan');
    mean_hit_fullBR = mean(hit_full_brs,'omitnan');
    mean_miss_fullBR = mean(miss_full_brs,'omitnan');
    
    allPLVs = [hitPLVs missPLVs];
    mean_all_PLVs = mean(allPLVs,'omitnan');
    mean_hit_PLVs = mean(hitPLVs,'omitnan');
    mean_miss_PLVs= mean(missPLVs,'omitnan');

    meanFastResp_ibi_s2 = mean(fastResp_ibi_s2,'omitnan'); meanSlowResp_ibi_s2 = mean(slowResp_ibi_s2,'omitnan');
    meanFastResp_ibi_s1 = mean(fastResp_ibi_s1,'omitnan'); meanSlowResp_ibi_s1 = mean(slowResp_ibi_s1,'omitnan');
    meanFastResp_ibi_s0 = mean(fastResp_ibi_s0,'omitnan'); meanSlowResp_ibi_s0 = mean(slowResp_ibi_s0,'omitnan');
    meanFastResp_ibi_s_p1 = mean(fastResp_ibi_s_p1,'omitnan'); meanSlowResp_ibi_s_p1 = mean(slowResp_ibi_s_p1,'omitnan');

    acd_inhale_s2 = acd_s2(inhale_s2);
    acd_inhaleHit_s2 = acd_s2(inhaleHit_s2);
    acd_inhaleMiss_s2 = acd_s2(inhaleMiss_s2);
    acd_inhale_s1 = acd_s1(inhale_s1);
    acd_inhaleHit_s1 = acd_s1(inhaleHit_s1);
    acd_inhaleMiss_s1 = acd_s1(inhaleMiss_s1);
    acd_inhale_s0 = acd_s0(inhale_s0);
    acd_inhaleHit_s0 = acd_s0(inhaleHit_s0);
    acd_inhaleMiss_s0 = acd_s0(inhaleMiss_s0);
    acd_inhale_s_p1 = acd_s_p1(inhale_s_p1);
    acd_inhaleHit_s_p1 = acd_s_p1(inhaleHit_s_p1);
    acd_inhaleMiss_s_p1 = acd_s_p1(inhaleMiss_s_p1);

    mean_acd_exhale_s2 = mean(acd_exhale_s2,'omitnan');
    mean_acd_exhaleHit_s2 = mean(acd_exhaleHit_s2,'omitnan');
    mean_acd_exhaleMiss_s2 = mean(acd_exhaleMiss_s2,'omitnan');
    mean_acd_exhale_s1 = mean(acd_exhale_s1,'omitnan');
    mean_acd_exhaleHit_s1 = mean(acd_exhaleHit_s1,'omitnan');
    mean_acd_exhaleMiss_s1 = mean(acd_exhaleMiss_s1,'omitnan');
    mean_acd_exhale_s0 = mean(acd_exhale_s0,'omitnan');
    mean_acd_exhaleHit_s0 = mean(acd_exhaleHit_s0,'omitnan');
    mean_acd_exhaleMiss_s0 = mean(acd_exhaleMiss_s0,'omitnan');
    mean_acd_exhale_s_p1 = mean(acd_exhale_s_p1,'omitnan');
    mean_acd_exhaleHit_s_p1 = mean(acd_exhaleHit_s_p1,'omitnan');
    mean_acd_exhaleMiss_s_p1 = mean(acd_exhaleMiss_s_p1,'omitnan');

    mean_acd_inhale_s2 = mean(acd_inhale_s2,'omitnan');
    mean_acd_inhaleHit_s2 = mean(acd_inhaleHit_s2,'omitnan');
    mean_acd_inhaleMiss_s2 = mean(acd_inhaleMiss_s2,'omitnan');
    mean_acd_inhale_s1 = mean(acd_inhale_s1,'omitnan');
    mean_acd_inhaleHit_s1 = mean(acd_inhaleHit_s1,'omitnan');
    mean_acd_inhaleMiss_s1 = mean(acd_inhaleMiss_s1,'omitnan');
    mean_acd_inhale_s0 = mean(acd_inhale_s0,'omitnan');
    mean_acd_inhaleHit_s0 = mean(acd_inhaleHit_s0,'omitnan');
    mean_acd_inhaleMiss_s0 = mean(acd_inhaleMiss_s0,'omitnan');
    mean_acd_inhale_s_p1 = mean(acd_inhale_s_p1,'omitnan');
    mean_acd_inhaleHit_s_p1 = mean(acd_inhaleHit_s_p1,'omitnan');
    mean_acd_inhaleMiss_s_p1 = mean(acd_inhaleMiss_s_p1,'omitnan');

    mean_exhalePeriod = mean(br_exhalePeriods,'omitnan');
    mean_inhalePeriod = mean(br_inhalePeriods,'omitnan');

    fileName_acd = 'acd_inhalationVsExhalation.mat'; 
    saveDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo];
    save(fullfile(saveDir, fileName_acd), 'mean_acd_exhale_s2', 'mean_acd_exhaleHit_s2', 'mean_acd_exhaleMiss_s2',...
        'mean_acd_exhale_s1', 'mean_acd_exhaleHit_s1', 'mean_acd_exhaleMiss_s1',...
        'mean_acd_exhale_s0', 'mean_acd_exhaleHit_s0', 'mean_acd_exhaleMiss_s0',...
        'mean_acd_exhale_s_p1', 'mean_acd_exhaleHit_s_p1', 'mean_acd_exhaleMiss_s_p1',...
        'mean_acd_inhale_s2', 'mean_acd_inhaleHit_s2', 'mean_acd_inhaleMiss_s2',...
        'mean_acd_inhale_s1', 'mean_acd_inhaleHit_s1', 'mean_acd_inhaleMiss_s1',...
        'mean_acd_inhale_s0', 'mean_acd_inhaleHit_s0', 'mean_acd_inhaleMiss_s0',...
        'mean_acd_inhale_s_p1', 'mean_acd_inhaleHit_s_p1', 'mean_acd_inhaleMiss_s_p1',...
        'mean_hitBrPeriods','mean_missBrPeriods','mean_hitExhaleBrPeriods','mean_missExhaleBrPeriods','mean_hitInhaleBrPeriods','mean_missInhaleBrPeriods',...
        'mean_all_brs','mean_hit_fullBR','mean_miss_fullBR',...
        'exh_s2_hitRate','inh_s2_hitRate','exh_s1_hitRate','inh_s1_hitRate','exh_s0_hitRate','inh_s0_hitRate','exh_stim_hitRate','inh_stim_hitRate',...
        'mean_all_PLVs','mean_hit_PLVs','mean_miss_PLVs');

    fileName_inh_exh_phases = 'acd_inhalationVsExhalationPhaseValues';
    saveDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo];
    save(fullfile(saveDir, fileName_inh_exh_phases),"inhale_s2Phase","inhale_s1Phase","inhale_s0Phase","exhale_s2Phase","exhale_s1Phase","exhale_s0Phase");

    %% plot circular to linear correlation plot for each subject!!!
    preStim_brRadians = deg2rad(preStim_brPhases);
    preStim_ACDs_val = preStim_ACDs(~isnan(preStim_ACDs) & ~isnan(preStim_brRadians));
    preStim_brRadians_val = preStim_brRadians(~isnan(preStim_ACDs) & ~isnan(preStim_brRadians));
    
    preStim_brRadians_hit = deg2rad(preStim_brPhases_hit);
    preStim_ACDs_hit_val = preStim_ACDs_hit(~isnan(preStim_ACDs_hit) & ~isnan(preStim_brRadians_hit));
    preStim_brRadians_hit_val = preStim_brRadians_hit(~isnan(preStim_ACDs_hit) & ~isnan(preStim_brRadians_hit));

    preStim_brRadians_miss = deg2rad(preStim_brPhases_miss);
    preStim_ACDs_miss_val = preStim_ACDs_miss(~isnan(preStim_ACDs_miss) & ~isnan(preStim_brRadians_miss));
    preStim_brRadians_miss_val = preStim_brRadians_miss(~isnan(preStim_ACDs_miss) & ~isnan(preStim_brRadians_miss));

    preStim_brRadians_fast = deg2rad(preStim_brPhases_fast);
    preStim_ACDs_fast_val = preStim_ACDs_fast(~isnan(preStim_ACDs_fast) & ~isnan(preStim_brRadians_fast));
    preStim_brRadians_fast_val = preStim_brRadians_fast(~isnan(preStim_ACDs_fast) & ~isnan(preStim_brRadians_fast));

    preStim_brRadians_slow = deg2rad(preStim_brPhases_slow);
    preStim_ACDs_slow_val = preStim_ACDs_slow(~isnan(preStim_ACDs_slow) & ~isnan(preStim_brRadians_slow));
    preStim_brRadians_slow_val = preStim_brRadians_slow(~isnan(preStim_ACDs_slow) & ~isnan(preStim_brRadians_slow));

    preStim_brRadians_highC = deg2rad(preStim_brPhases_highC);
    preStim_ACDs_highC_val = preStim_ACDs_highC(~isnan(preStim_ACDs_highC) & ~isnan(preStim_brRadians_highC));
    preStim_brRadians_highC_val = preStim_brRadians_highC(~isnan(preStim_ACDs_highC) & ~isnan(preStim_brRadians_highC));

    preStim_brRadians_lowC = deg2rad(preStim_brPhases_lowC);
    preStim_ACDs_lowC_val = preStim_ACDs_lowC(~isnan(preStim_ACDs_lowC) & ~isnan(preStim_brRadians_lowC));
    preStim_brRadians_lowC_val = preStim_brRadians_lowC(~isnan(preStim_ACDs_lowC) & ~isnan(preStim_brRadians_lowC));

    fastACD = mean(preStim_ACDs_fast_val,'omitnan');
    slowACD = mean(preStim_ACDs_slow_val,'omitnan');

    %% save the breathing phases and delta IBIs corresponding to these phases!
    fileName_circLinear = 'acd_preStim_brPhase.mat';
    saveDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo];
    save(fullfile(saveDir,fileName_circLinear),'preStim_brRadians_val','preStim_brRadians_hit_val','preStim_brRadians_miss_val','preStim_brRadians_highC_val','preStim_brRadians_lowC_val',...
        'preStim_ACDs_val','preStim_ACDs_hit_val','preStim_ACDs_miss_val','preStim_ACDs_highC_val','preStim_ACDs_lowC_val',...
        'preStim_brRadians_fast_val','preStim_brRadians_slow_val',...
        'preStim_ACDs_fast_val','preStim_ACDs_slow_val',...
        'acd_s2','acd_s1','acd_s0','acd_s_p1',...
        'acd_s2_hit','acd_s1_hit','acd_s0_hit','acd_s_p1_hit','acd_s2_highC','acd_s1_highC','acd_s0_highC',...
        'acd_s2_miss','acd_s1_miss','acd_s0_miss','acd_s_p1_miss','acd_s2_lowC','acd_s1_lowC','acd_s0_lowC');

    fileName_fastSlowResp_ibi = 'ibi_fastSlowResp.mat';
    save(fullfile(saveDir,fileName_fastSlowResp_ibi),'meanFastResp_ibi_s2','meanSlowResp_ibi_s2','meanFastResp_ibi_s1','meanSlowResp_ibi_s1',...
        'meanFastResp_ibi_s0','meanSlowResp_ibi_s0','meanFastResp_ibi_s_p1','meanSlowResp_ibi_s_p1');

    %% fit sinusoidal model to the brPhase ~ deltaIBI polar distribution
    addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\code\VTD_rsaSinusoidalModel
    % Example data (replace with your actual vectors)
    % phase_rad = ...     % Nx1 vector of breathing phase (radians)
    % amplitude = ...     % Nx1 vector of response magnitude (e.g., ΔIBI)

    [A, B, phi, y, SSres, SStot, Rsq, phase, amp, fitX, fitY] = RSA_sinusoidalFit(preStim_brRadians_val',preStim_ACDs_val');

    %% fit only for s-2 / s-1 / s-0 heartbeats
    %%% s3 all
    s3_brRadians = deg2rad(s_min3_allPhases);
    s3_ACDs_val = acd_s3(~isnan(acd_s3') & ~isnan(s3_brRadians));
    s3_brRadians_val = s3_brRadians(~isnan(acd_s3') & ~isnan(s3_brRadians));
    %%% s3 hit
    s3_brRadians_hit = deg2rad(s_min3_hitPhases);
    s3_ACDs_hit_val = acd_s3_hit(~isnan(acd_s3_hit) & ~isnan(s3_brRadians_hit));
    s3_brRadians_hit_val = s3_brRadians_hit(~isnan(acd_s3_hit) & ~isnan(s3_brRadians_hit));
    %%% s3 miss
    s3_brRadians_miss = deg2rad(s_min3_missPhases);
    s3_ACDs_miss_val = acd_s3_miss(~isnan(acd_s3_miss) & ~isnan(s3_brRadians_miss));
    s3_brRadians_miss_val = s3_brRadians_miss(~isnan(acd_s3_miss) & ~isnan(s3_brRadians_miss));
    
    %%% s2 all
    s2_brRadians = deg2rad(s_min2_allPhases);
    s2_ACDs_val = acd_s2(~isnan(acd_s2') & ~isnan(s2_brRadians));
    s2_brRadians_val = s2_brRadians(~isnan(acd_s2') & ~isnan(s2_brRadians));
    %%% s2 hit
    s2_brRadians_hit = deg2rad(s_min2_hitPhases);
    s2_ACDs_hit_val = acd_s2_hit(~isnan(acd_s2_hit) & ~isnan(s2_brRadians_hit));
    s2_brRadians_hit_val = s2_brRadians_hit(~isnan(acd_s2_hit) & ~isnan(s2_brRadians_hit));
    %%% s2 miss
    s2_brRadians_miss = deg2rad(s_min2_missPhases);
    s2_ACDs_miss_val = acd_s2_miss(~isnan(acd_s2_miss) & ~isnan(s2_brRadians_miss));
    s2_brRadians_miss_val = s2_brRadians_miss(~isnan(acd_s2_miss) & ~isnan(s2_brRadians_miss));
    %%% s2 highC
    s2_brRadians_highC = deg2rad(s_min2_highCPhases);
    s2_ACDs_highC_val = acd_s2_highC(~isnan(acd_s2_highC) & ~isnan(s2_brRadians_highC));
    s2_brRadians_highC_val = s2_brRadians_highC(~isnan(acd_s2_highC) & ~isnan(s2_brRadians_highC));
    %%% s2 lowC
    s2_brRadians_lowC = deg2rad(s_min2_lowCPhases);
    s2_ACDs_lowC_val = acd_s2_lowC(~isnan(acd_s2_lowC) & ~isnan(s2_brRadians_lowC));
    s2_brRadians_lowC_val = s2_brRadians_lowC(~isnan(acd_s2_lowC) & ~isnan(s2_brRadians_lowC));
    %%% s2 HH
    s2_brRadians_HH = deg2rad(s_min2_HHPhases);
    s2_ACDs_HH_val = acd_s2_HH(~isnan(acd_s2_HH) & ~isnan(s2_brRadians_HH));
    s2_brRadians_HH_val = s2_brRadians_HH(~isnan(acd_s2_HH) & ~isnan(s2_brRadians_HH));
    %%% s2 HL
    s2_brRadians_HL = deg2rad(s_min2_HLPhases);
    s2_ACDs_HL_val = acd_s2_HL(~isnan(acd_s2_HL) & ~isnan(s2_brRadians_HL));
    s2_brRadians_HL_val = s2_brRadians_HL(~isnan(acd_s2_HL) & ~isnan(s2_brRadians_HL));
    %%% s2 MH
    s2_brRadians_MH = deg2rad(s_min2_MHPhases);
    s2_ACDs_MH_val = acd_s2_MH(~isnan(acd_s2_MH) & ~isnan(s2_brRadians_MH));
    s2_brRadians_MH_val = s2_brRadians_MH(~isnan(acd_s2_MH) & ~isnan(s2_brRadians_MH));
    %%% s2 ML
    s2_brRadians_ML = deg2rad(s_min2_MLPhases);
    s2_ACDs_ML_val = acd_s2_ML(~isnan(acd_s2_ML) & ~isnan(s2_brRadians_ML));
    s2_brRadians_ML_val = s2_brRadians_ML(~isnan(acd_s2_ML) & ~isnan(s2_brRadians_ML));

    %%% s1 all
    s1_brRadians = deg2rad(s_min1_allPhases);
    s1_ACDs_val = acd_s1(~isnan(acd_s1') & ~isnan(s1_brRadians));
    s1_brRadians_val = s1_brRadians(~isnan(acd_s1') & ~isnan(s1_brRadians));
    %%% s1 hit
    s1_brRadians_hit = deg2rad(s_min1_hitPhases);
    s1_ACDs_hit_val = acd_s1_hit(~isnan(acd_s1_hit) & ~isnan(s1_brRadians_hit));
    s1_brRadians_hit_val = s1_brRadians_hit(~isnan(acd_s1_hit) & ~isnan(s1_brRadians_hit));
    %%% s1 miss
    s1_brRadians_miss = deg2rad(s_min1_missPhases);
    s1_ACDs_miss_val = acd_s1_miss(~isnan(acd_s1_miss) & ~isnan(s1_brRadians_miss));
    s1_brRadians_miss_val = s1_brRadians_miss(~isnan(acd_s1_miss) & ~isnan(s1_brRadians_miss));
    %%% s1 highC
    s1_brRadians_highC = deg2rad(s_min1_highCPhases);
    s1_ACDs_highC_val = acd_s1_highC(~isnan(acd_s1_highC) & ~isnan(s1_brRadians_highC));
    s1_brRadians_highC_val = s1_brRadians_highC(~isnan(acd_s1_highC) & ~isnan(s1_brRadians_highC));
    %%% s1 lowC
    s1_brRadians_lowC = deg2rad(s_min1_lowCPhases);
    s1_ACDs_lowC_val = acd_s1_lowC(~isnan(acd_s1_lowC) & ~isnan(s1_brRadians_lowC));
    s1_brRadians_lowC_val = s1_brRadians_lowC(~isnan(acd_s1_lowC) & ~isnan(s1_brRadians_lowC));

    %%% s0 all
    s0_brRadians = deg2rad(s_minZ_allPhases);
    s0_ACDs_val = acd_s0(~isnan(acd_s0') & ~isnan(s0_brRadians));
    s0_brRadians_val = s0_brRadians(~isnan(acd_s0') & ~isnan(s0_brRadians));
    %%% s0 hit
    s0_brRadians_hit = deg2rad(s_minZ_hitPhases);
    s0_ACDs_hit_val = acd_s0_hit(~isnan(acd_s0_hit) & ~isnan(s0_brRadians_hit));
    s0_brRadians_hit_val = s0_brRadians_hit(~isnan(acd_s0_hit) & ~isnan(s0_brRadians_hit));
    %%% s0 miss
    s0_brRadians_miss = deg2rad(s_minZ_missPhases);
    s0_ACDs_miss_val = acd_s0_miss(~isnan(acd_s0_miss) & ~isnan(s0_brRadians_miss));
    s0_brRadians_miss_val = s0_brRadians_miss(~isnan(acd_s0_miss) & ~isnan(s0_brRadians_miss));
    %%% s0 highC
    s0_brRadians_highC = deg2rad(s_minZ_highCPhases);
    s0_ACDs_highC_val = acd_s0_highC(~isnan(acd_s0_highC) & ~isnan(s0_brRadians_highC));
    s0_brRadians_highC_val = s0_brRadians_highC(~isnan(acd_s0_highC) & ~isnan(s0_brRadians_highC));
    %%% s0 lowC
    s0_brRadians_lowC = deg2rad(s_minZ_lowCPhases);
    s0_ACDs_lowC_val = acd_s0_lowC(~isnan(acd_s0_lowC) & ~isnan(s0_brRadians_lowC));
    s0_brRadians_lowC_val = s0_brRadians_lowC(~isnan(acd_s0_lowC) & ~isnan(s0_brRadians_lowC));
    %%% s0 HH
    s0_brRadians_HH = deg2rad(s_minZ_HHPhases);
    s0_ACDs_HH_val = acd_s0_HH(~isnan(acd_s0_HH) & ~isnan(s0_brRadians_HH));
    s0_brRadians_HH_val = s0_brRadians_HH(~isnan(acd_s0_HH) & ~isnan(s0_brRadians_HH));
    %%% s0 HL
    s0_brRadians_HL = deg2rad(s_minZ_HLPhases);
    s0_ACDs_HL_val = acd_s0_HL(~isnan(acd_s0_HL) & ~isnan(s0_brRadians_HL));
    s0_brRadians_HL_val = s0_brRadians_HL(~isnan(acd_s0_HL) & ~isnan(s0_brRadians_HL));
    %%% s0 MH
    s0_brRadians_MH = deg2rad(s_minZ_MHPhases);
    s0_ACDs_MH_val = acd_s0_MH(~isnan(acd_s0_MH) & ~isnan(s0_brRadians_MH));
    s0_brRadians_MH_val = s0_brRadians_MH(~isnan(acd_s0_MH) & ~isnan(s0_brRadians_MH));
    %%% s0 ML
    s0_brRadians_ML = deg2rad(s_minZ_MLPhases);
    s0_ACDs_ML_val = acd_s0_ML(~isnan(acd_s0_ML) & ~isnan(s0_brRadians_ML));
    s0_brRadians_ML_val = s0_brRadians_ML(~isnan(acd_s0_ML) & ~isnan(s0_brRadians_ML));

    %%% s_p1 all
    s_p1_brRadians = deg2rad(s_plus1_allPhases);
    s_p1_ACDs_val = acd_s_p1(~isnan(acd_s_p1') & ~isnan(s_p1_brRadians));
    s_p1_brRadians_val = s_p1_brRadians(~isnan(acd_s_p1') & ~isnan(s_p1_brRadians));
    %%% s_p1 hit
    s_p1_brRadians_hit = deg2rad(s_p1_hitPhases);
    s_p1_ACDs_hit_val = acd_s_p1_hit(~isnan(acd_s_p1_hit) & ~isnan(s_p1_brRadians_hit));
    s_p1_brRadians_hit_val = s_p1_brRadians_hit(~isnan(acd_s_p1_hit) & ~isnan(s_p1_brRadians_hit));
    %%% s_p1 miss
    s_p1_brRadians_miss = deg2rad(s_p1_missPhases);
    s_p1_ACDs_miss_val = acd_s_p1_miss(~isnan(acd_s_p1_miss) & ~isnan(s_p1_brRadians_miss));
    s_p1_brRadians_miss_val = s_p1_brRadians_miss(~isnan(acd_s_p1_miss) & ~isnan(s_p1_brRadians_miss));

    %% separate fits beat-by-beat
    %%% s-3 fit
    [A_s3, B_s3, phi_s3, y_s3, SSres_s3, SStot_s3, Rsq_s3, phase_s3, amp_s3, fitX_s3, fitY_s3] = RSA_sinusoidalFit(s3_brRadians_val,s3_ACDs_val');

    %%% s-2 fit
    [A_s2, B_s2, phi_s2, y_s2, SSres_s2, SStot_s2, Rsq_s2, phase_s2, amp_s2, fitX_s2, fitY_s2] = RSA_sinusoidalFit(s2_brRadians_val,s2_ACDs_val');

    %%% s-1 fit
    [A_s1, B_s1, phi_s1, y_s1, SSres_s1, SStot_s1, Rsq_s1, phase_s1, amp_s1, fitX_s1, fitY_s1] = RSA_sinusoidalFit(s1_brRadians_val,s1_ACDs_val');

    %%% s-0 fit
    [A_s0, B_s0, phi_s0, y_s0, SSres_s0, SStot_s0, Rsq_s0, phase_s0, amp_s0, fitX_s0, fitY_s0] = RSA_sinusoidalFit(s0_brRadians_val,s0_ACDs_val');

    %%% s-p1 fit
    [A_s_p1, B_s_p1, phi_s_p1, y_s_p1, SSres_s_p1, SStot_s_p1, Rsq_s_p1, phase_s_p1, amp_s_p1, fitX_s_p1, fitY_s_p1] = RSA_sinusoidalFit(s_p1_brRadians_val,s_p1_ACDs_val');

    %% separate fits for hit and miss trials (avg pre stim)
    %%% ONLY HITS
    [A_hit, B_hit, phi_hit, y_hit, SSres_hit, SStot_hit, Rsq_hit, phase_hit, amp_hit, fitX_hit, fitY_hit] = RSA_sinusoidalFit(preStim_brRadians_hit_val',preStim_ACDs_hit_val');

    %%% ONLY MISSES
    [A_miss, B_miss, phi_miss, y_miss, SSres_miss, SStot_miss, Rsq_miss, phase_miss, amp_miss, fitX_miss, fitY_miss] = RSA_sinusoidalFit(preStim_brRadians_miss_val',preStim_ACDs_miss_val');
    if ~ismember(subjectNo,{'05','16','26','28'})
        %%% ONLY HIGH-C
        [A_highC, B_highC, phi_highC, y_highC, SSres_highC, SStot_highC, Rsq_highC, phase_highC, amp_highC, fitX_highC, fitY_highC] = RSA_sinusoidalFit(preStim_brRadians_highC_val',preStim_ACDs_highC_val');
    
        %%% ONLY LOW-C
        [A_lowC, B_lowC, phi_lowC, y_lowC, SSres_lowC, SStot_lowC, Rsq_lowC, phase_lowC, amp_lowC, fitX_lowC, fitY_lowC] = RSA_sinusoidalFit(preStim_brRadians_lowC_val',preStim_ACDs_lowC_val');
    else
       fitY_highC=NaN; fitX_highC= fitY_highC; amp_highC= fitX_highC; phase_highC= amp_highC; Rsq_highC= phase_highC; SStot_highC= Rsq_highC; SSres_highC= SStot_highC; y_highC= SSres_highC; phi_highC= y_highC; B_highC= phi_highC; A_highC= B_highC;
       fitY_lowC=NaN; fitX_lowC= fitY_lowC; amp_lowC= fitX_lowC; phase_lowC= amp_lowC; Rsq_lowC= phase_lowC; SStot_lowC= Rsq_lowC; SSres_lowC= SStot_lowC; y_lowC= SSres_lowC; phi_lowC= y_lowC; B_lowC= phi_lowC; A_lowC= B_lowC;
    end
    %% separate fits for hit and miss trials: beat-by-beat
    %%% HITS
    [A_hit_s3, B_hit_s3, phi_hit_s3, y_hit_s3, SSres_hit_s3, SStot_hit_s3, Rsq_hit_s3, phase_hit_s3, amp_hit_s3, fitX_hit_s3, fitY_hit_s3] = RSA_sinusoidalFit(s3_brRadians_hit_val',s3_ACDs_hit_val');

    [A_hit_s2, B_hit_s2, phi_hit_s2, y_hit_s2, SSres_hit_s2, SStot_hit_s2, Rsq_hit_s2, phase_hit_s2, amp_hit_s2, fitX_hit_s2, fitY_hit_s2] = RSA_sinusoidalFit(s2_brRadians_hit_val',s2_ACDs_hit_val');

    [A_hit_s1, B_hit_s1, phi_hit_s1, y_hit_s1, SSres_hit_s1, SStot_hit_s1, Rsq_hit_s1, phase_hit_s1, amp_hit_s1, fitX_hit_s1, fitY_hit_s1] = RSA_sinusoidalFit(s1_brRadians_hit_val',s1_ACDs_hit_val');

    [A_hit_s0, B_hit_s0, phi_hit_s0, y_hit_s0, SSres_hit_s0, SStot_hit_s0, Rsq_hit_s0, phase_hit_s0, amp_hit_s0, fitX_hit_s0, fitY_hit_s0] = RSA_sinusoidalFit(s0_brRadians_hit_val',s0_ACDs_hit_val');

    [A_hit_s_p1, B_hit_s_p1, phi_hit_s_p1, y_hit_s_p1, SSres_hit_s_p1, SStot_hit_s_p1, Rsq_hit_s_p1, phase_hit_s_p1, amp_hit_s_p1, fitX_hit_s_p1, fitY_hit_s_p1] = RSA_sinusoidalFit(s_p1_brRadians_hit_val',s_p1_ACDs_hit_val');

    %%% MISSES
    [A_miss_s3, B_miss_s3, phi_miss_s3, y_miss_s3, SSres_miss_s3, SStot_miss_s3, Rsq_miss_s3, phase_miss_s3, amp_miss_s3, fitX_miss_s3, fitY_miss_s3] = RSA_sinusoidalFit(s3_brRadians_miss_val',s3_ACDs_miss_val');
    
    [A_miss_s2, B_miss_s2, phi_miss_s2, y_miss_s2, SSres_miss_s2, SStot_miss_s2, Rsq_miss_s2, phase_miss_s2, amp_miss_s2, fitX_miss_s2, fitY_miss_s2] = RSA_sinusoidalFit(s2_brRadians_miss_val',s2_ACDs_miss_val');

    [A_miss_s1, B_miss_s1, phi_miss_s1, y_miss_s1, SSres_miss_s1, SStot_miss_s1, Rsq_miss_s1, phase_miss_s1, amp_miss_s1, fitX_miss_s1, fitY_miss_s1] = RSA_sinusoidalFit(s1_brRadians_miss_val',s1_ACDs_miss_val');

    [A_miss_s0, B_miss_s0, phi_miss_s0, y_miss_s0, SSres_miss_s0, SStot_miss_s0, Rsq_miss_s0, phase_miss_s0, amp_miss_s0, fitX_miss_s0, fitY_miss_s0] = RSA_sinusoidalFit(s0_brRadians_miss_val',s0_ACDs_miss_val');

    [A_miss_s_p1, B_miss_s_p1, phi_miss_s_p1, y_miss_s_p1, SSres_miss_s_p1, SStot_miss_s_p1, Rsq_miss_s_p1, phase_miss_s_p1, amp_miss_s_p1, fitX_miss_s_p1, fitY_miss_s_p1] = RSA_sinusoidalFit(s_p1_brRadians_miss_val',s_p1_ACDs_miss_val');
    if ~ismember(subjectNo,{'05','16','26','28'})
        %%% HIGH-C
        [A_highC_s2, B_highC_s2, phi_highC_s2, y_highC_s2, SSres_highC_s2, SStot_highC_s2, Rsq_highC_s2, phase_highC_s2, amp_highC_s2, fitX_highC_s2, fitY_highC_s2] = RSA_sinusoidalFit(s2_brRadians_highC_val',s2_ACDs_highC_val');
    
        [A_highC_s1, B_highC_s1, phi_highC_s1, y_highC_s1, SSres_highC_s1, SStot_highC_s1, Rsq_highC_s1, phase_highC_s1, amp_highC_s1, fitX_highC_s1, fitY_highC_s1] = RSA_sinusoidalFit(s1_brRadians_highC_val',s1_ACDs_highC_val');
    
        [A_highC_s0, B_highC_s0, phi_highC_s0, y_highC_s0, SSres_highC_s0, SStot_highC_s0, Rsq_highC_s0, phase_highC_s0, amp_highC_s0, fitX_highC_s0, fitY_highC_s0] = RSA_sinusoidalFit(s0_brRadians_highC_val',s0_ACDs_highC_val');
    
        %%% LOW-C
        [A_lowC_s2, B_lowC_s2, phi_lowC_s2, y_lowC_s2, SSres_lowC_s2, SStot_lowC_s2, Rsq_lowC_s2, phase_lowC_s2, amp_lowC_s2, fitX_lowC_s2, fitY_lowC_s2] = RSA_sinusoidalFit(s2_brRadians_lowC_val',s2_ACDs_lowC_val');
    
        [A_lowC_s1, B_lowC_s1, phi_lowC_s1, y_lowC_s1, SSres_lowC_s1, SStot_lowC_s1, Rsq_lowC_s1, phase_lowC_s1, amp_lowC_s1, fitX_lowC_s1, fitY_lowC_s1] = RSA_sinusoidalFit(s1_brRadians_lowC_val',s1_ACDs_lowC_val');
    
        [A_lowC_s0, B_lowC_s0, phi_lowC_s0, y_lowC_s0, SSres_lowC_s0, SStot_lowC_s0, Rsq_lowC_s0, phase_lowC_s0, amp_lowC_s0, fitX_lowC_s0, fitY_lowC_s0] = RSA_sinusoidalFit(s0_brRadians_lowC_val',s0_ACDs_lowC_val');
    else
        fitY_highC_s2=NaN; fitX_highC_s2= fitY_highC_s2; amp_highC_s2= fitX_highC_s2; phase_highC_s2= amp_highC_s2; Rsq_highC_s2= phase_highC_s2; SStot_highC_s2= Rsq_highC_s2; SSres_highC_s2= SStot_highC_s2; y_highC_s2= SSres_highC_s2; phi_highC_s2= y_highC_s2; B_highC_s2= phi_highC_s2; A_highC_s2= B_highC_s2;
        fitY_lowC_s2=NaN; fitX_lowC_s2= fitY_lowC_s2; amp_lowC_s2= fitX_lowC_s2; phase_lowC_s2= amp_lowC_s2; Rsq_lowC_s2= phase_lowC_s2; SStot_lowC_s2= Rsq_lowC_s2; SSres_lowC_s2= SStot_lowC_s2; y_lowC_s2= SSres_lowC_s2; phi_lowC_s2= y_lowC_s2; B_lowC_s2= phi_lowC_s2; A_lowC_s2= B_lowC_s2;
        fitY_highC_s1=NaN; fitX_highC_s1= fitY_highC_s1; amp_highC_s1= fitX_highC_s1; phase_highC_s1= amp_highC_s1; Rsq_highC_s1= phase_highC_s1; SStot_highC_s1= Rsq_highC_s1; SSres_highC_s1= SStot_highC_s1; y_highC_s1= SSres_highC_s1; phi_highC_s1= y_highC_s1; B_highC_s1= phi_highC_s1; A_highC_s1= B_highC_s1;
        fitY_lowC_s1=NaN; fitX_lowC_s1= fitY_lowC_s1; amp_lowC_s1= fitX_lowC_s1; phase_lowC_s1= amp_lowC_s1; Rsq_lowC_s1= phase_lowC_s1; SStot_lowC_s1= Rsq_lowC_s1; SSres_lowC_s1= SStot_lowC_s1; y_lowC_s1= SSres_lowC_s1; phi_lowC_s1= y_lowC_s1; B_lowC_s1= phi_lowC_s1; A_lowC_s1= B_lowC_s1;
        fitY_highC_s0=NaN; fitX_highC_s0= fitY_highC_s0; amp_highC_s0= fitX_highC_s0; phase_highC_s0= amp_highC_s0; Rsq_highC_s0= phase_highC_s0; SStot_highC_s0= Rsq_highC_s0; SSres_highC_s0= SStot_highC_s0; y_highC_s0= SSres_highC_s0; phi_highC_s0= y_highC_s0; B_highC_s0= phi_highC_s0; A_highC_s0= B_highC_s0;
        fitY_lowC_s0=NaN; fitX_lowC_s0= fitY_lowC_s0; amp_lowC_s0= fitX_lowC_s0; phase_lowC_s0= amp_lowC_s0; Rsq_lowC_s0= phase_lowC_s0; SStot_lowC_s0= Rsq_lowC_s0; SSres_lowC_s0= SStot_lowC_s0; y_lowC_s0= SSres_lowC_s0; phi_lowC_s0= y_lowC_s0; B_lowC_s0= phi_lowC_s0; A_lowC_s0= B_lowC_s0;
    end
    %%% HH
%     [A_HH_s2, B_HH_s2, phi_HH_s2, y_HH_s2, SSres_HH_s2, SStot_HH_s2, Rsq_HH_s2, phase_HH_s2, amp_HH_s2, fitX_HH_s2, fitY_HH_s2] = RSA_sinusoidalFit(s2_brRadians_HH_val',s2_ACDs_HH_val');
%     
%     [A_HH_s0, B_HH_s0, phi_HH_s0, y_HH_s0, SSres_HH_s0, SStot_HH_s0, Rsq_HH_s0, phase_HH_s0, amp_HH_s0, fitX_HH_s0, fitY_HH_s0] = RSA_sinusoidalFit(s0_brRadians_HH_val',s0_ACDs_HH_val');
% 
%     %%% HL
%     [A_HL_s2, B_HL_s2, phi_HL_s2, y_HL_s2, SSres_HL_s2, SStot_HL_s2, Rsq_HL_s2, phase_HL_s2, amp_HL_s2, fitX_HL_s2, fitY_HL_s2] = RSA_sinusoidalFit(s2_brRadians_HL_val',s2_ACDs_HL_val');
%     
%     [A_HL_s0, B_HL_s0, phi_HL_s0, y_HL_s0, SSres_HL_s0, SStot_HL_s0, Rsq_HL_s0, phase_HL_s0, amp_HL_s0, fitX_HL_s0, fitY_HL_s0] = RSA_sinusoidalFit(s0_brRadians_HL_val',s0_ACDs_HL_val');
% 
%     %%% MH
%     [A_MH_s2, B_MH_s2, phi_MH_s2, y_MH_s2, SSres_MH_s2, SStot_MH_s2, Rsq_MH_s2, phase_MH_s2, amp_MH_s2, fitX_MH_s2, fitY_MH_s2] = RSA_sinusoidalFit(s2_brRadians_MH_val',s2_ACDs_MH_val');
%     
%     [A_MH_s0, B_MH_s0, phi_MH_s0, y_MH_s0, SSres_MH_s0, SStot_MH_s0, Rsq_MH_s0, phase_MH_s0, amp_MH_s0, fitX_MH_s0, fitY_MH_s0] = RSA_sinusoidalFit(s0_brRadians_MH_val',s0_ACDs_MH_val');
% 
%     %%% ML
%     [A_ML_s2, B_ML_s2, phi_ML_s2, y_ML_s2, SSres_ML_s2, SStot_ML_s2, Rsq_ML_s2, phase_ML_s2, amp_ML_s2, fitX_ML_s2, fitY_ML_s2] = RSA_sinusoidalFit(s2_brRadians_ML_val',s2_ACDs_ML_val');
%     
%     [A_ML_s0, B_ML_s0, phi_ML_s0, y_ML_s0, SSres_ML_s0, SStot_ML_s0, Rsq_ML_s0, phase_ML_s0, amp_ML_s0, fitX_ML_s0, fitY_ML_s0] = RSA_sinusoidalFit(s0_brRadians_ML_val',s0_ACDs_ML_val');
    

    %% separate fits for fast and slow response trials (avg pre stim)
    %%% ONLY FAST
    [A_fast, B_fast, phi_fast, y_fast, SSres_fast, SStot_fast, Rsq_fast, phase_fast, amp_fast, fitX_fast, fitY_fast] = RSA_sinusoidalFit(preStim_brRadians_fast_val',preStim_ACDs_fast_val');

    %%% ONLY SLOWS
    [A_slow, B_slow, phi_slow, y_slow, SSres_slow, SStot_slow, Rsq_slow, phase_slow, amp_slow, fitX_slow, fitY_slow] = RSA_sinusoidalFit(preStim_brRadians_slow_val',preStim_ACDs_slow_val');

    %% save sinusoidal model parameters
    filename_sineModel = 'rsa_ibi_sineModel.mat';
    save(fullfile(saveDir, filename_sineModel), 'A','A_hit','A_miss','A_highC','A_lowC','A_fast','A_slow', 'phi','phi_hit','phi_miss','phi_highC','phi_lowC','phi_fast','phi_slow',...
        'B','B_hit','B_miss','B_highC','B_lowC','B_fast','B_slow','y','y_hit','y_miss','y_highC','y_lowC','y_fast','y_slow','SSres','SSres_hit','SSres_miss','SSres_highC','SSres_lowC','SSres_fast','SSres_slow',...
        'SStot','SStot_hit','SStot_miss','SStot_highC','SStot_lowC','SStot_fast','SStot_slow','Rsq','Rsq_hit','Rsq_miss','Rsq_highC','Rsq_lowC','Rsq_fast','Rsq_slow','fastACD','slowACD');

    filename_sineEachBeat = 'rsa_ibi_sineModel_beatByBeat.mat';
    save(fullfile(saveDir,filename_sineEachBeat), 'A_s3', 'B_s3', 'phi_s3', 'y_s3', 'SSres_s3', 'SStot_s3', 'Rsq_s3', 'phase_s3', 'amp_s3', 'fitX_s3', 'fitY_s3',...
        'A_s2', 'B_s2', 'phi_s2', 'y_s2', 'SSres_s2', 'SStot_s2', 'Rsq_s2', 'phase_s2', 'amp_s2', 'fitX_s2', 'fitY_s2',...
        'A_s1', 'B_s1', 'phi_s1', 'y_s1', 'SSres_s1', 'SStot_s1', 'Rsq_s1', 'phase_s1', 'amp_s1', 'fitX_s1', 'fitY_s1',...
        'A_s0', 'B_s0', 'phi_s0', 'y_s0', 'SSres_s0', 'SStot_s0', 'Rsq_s0', 'phase_s0', 'amp_s0', 'fitX_s0', 'fitY_s0',...
        'A_s_p1', 'B_s_p1', 'phi_s_p1', 'y_s_p1', 'SSres_s_p1', 'SStot_s_p1', 'Rsq_s_p1', 'phase_s_p1', 'amp_s_p1', 'fitX_s_p1', 'fitY_s_p1',...
        'A_hit_s3','B_hit_s3', 'phi_hit_s3', 'y_hit_s3', 'SSres_hit_s3', 'SStot_hit_s3', 'Rsq_hit_s3','fitX_hit_s3','fitY_hit_s3',...
        'A_hit_s2','B_hit_s2', 'phi_hit_s2', 'y_hit_s2', 'SSres_hit_s2', 'SStot_hit_s2', 'Rsq_hit_s2','fitX_hit_s2','fitY_hit_s2',...
        'A_hit_s1','B_hit_s1', 'phi_hit_s1', 'y_hit_s1', 'SSres_hit_s1', 'SStot_hit_s1', 'Rsq_hit_s1','fitX_hit_s1','fitY_hit_s1',...
        'A_hit_s0','B_hit_s0', 'phi_hit_s0', 'y_hit_s0', 'SSres_hit_s0', 'SStot_hit_s0', 'Rsq_hit_s0','fitX_hit_s0','fitY_hit_s0',...
        'A_hit_s_p1','B_hit_s_p1', 'phi_hit_s_p1', 'y_hit_s_p1', 'SSres_hit_s_p1', 'SStot_hit_s_p1', 'Rsq_hit_s_p1',...
        'A_miss_s3','B_miss_s3', 'phi_miss_s3', 'y_miss_s3', 'SSres_miss_s3', 'SStot_miss_s3', 'Rsq_miss_s3','fitX_miss_s3','fitY_miss_s3',...
        'A_miss_s2','B_miss_s2', 'phi_miss_s2', 'y_miss_s2', 'SSres_miss_s2', 'SStot_miss_s2', 'Rsq_miss_s2','fitX_miss_s2','fitY_miss_s2',...
        'A_miss_s1','B_miss_s1', 'phi_miss_s1', 'y_miss_s1', 'SSres_miss_s1', 'SStot_miss_s1', 'Rsq_miss_s1','fitX_miss_s1','fitY_miss_s1',...
        'A_miss_s0','B_miss_s0', 'phi_miss_s0', 'y_miss_s0', 'SSres_miss_s0', 'SStot_miss_s0', 'Rsq_miss_s0','fitX_miss_s0','fitY_miss_s0',...
        'A_miss_s_p1','B_miss_s_p1', 'phi_miss_s_p1', 'y_miss_s_p1', 'SSres_miss_s_p1', 'SStot_miss_s_p1', 'Rsq_miss_s_p1',...
        'A_highC_s2','B_highC_s2', 'phi_highC_s2', 'y_highC_s2', 'SSres_highC_s2', 'SStot_highC_s2', 'Rsq_highC_s2','fitX_highC_s2','fitY_highC_s2',...
        'A_highC_s1','B_highC_s1', 'phi_highC_s1', 'y_highC_s1', 'SSres_highC_s1', 'SStot_highC_s1', 'Rsq_highC_s1','fitX_highC_s1','fitY_highC_s1',...
        'A_highC_s0','B_highC_s0', 'phi_highC_s0', 'y_highC_s0', 'SSres_highC_s0', 'SStot_highC_s0', 'Rsq_highC_s0','fitX_highC_s0','fitY_highC_s0',...
        'A_lowC_s2','B_lowC_s2', 'phi_lowC_s2', 'y_lowC_s2', 'SSres_lowC_s2', 'SStot_lowC_s2', 'Rsq_lowC_s2','fitX_lowC_s2','fitY_lowC_s2',...
        'A_lowC_s1','B_lowC_s1', 'phi_lowC_s1', 'y_lowC_s1', 'SSres_lowC_s1', 'SStot_lowC_s1', 'Rsq_lowC_s1','fitX_lowC_s1','fitY_lowC_s1',...
        'A_lowC_s0','B_lowC_s0', 'phi_lowC_s0', 'y_lowC_s0', 'SSres_lowC_s0', 'SStot_lowC_s0', 'Rsq_lowC_s0','fitX_lowC_s0','fitY_lowC_s0');

%     filename_2x2 = 'rsa_ibi_sineModel_2x2.mat';
%     save(fullfile(saveDir,filename_2x2),'A_HH_s2', 'B_HH_s2', 'phi_HH_s2', 'y_HH_s2', 'SSres_HH_s2', 'SStot_HH_s2', 'Rsq_HH_s2', 'phase_HH_s2', 'amp_HH_s2', 'fitX_HH_s2', 'fitY_HH_s2',...
%         'A_HL_s2', 'B_HL_s2', 'phi_HL_s2', 'y_HL_s2', 'SSres_HL_s2', 'SStot_HL_s2', 'Rsq_HL_s2', 'phase_HL_s2', 'amp_HL_s2', 'fitX_HL_s2', 'fitY_HL_s2',...
%         'A_MH_s2', 'B_MH_s2', 'phi_MH_s2', 'y_MH_s2', 'SSres_MH_s2', 'SStot_MH_s2', 'Rsq_MH_s2', 'phase_MH_s2', 'amp_MH_s2', 'fitX_MH_s2', 'fitY_MH_s2',...
%         'A_ML_s2', 'B_ML_s2', 'phi_ML_s2', 'y_ML_s2', 'SSres_ML_s2', 'SStot_ML_s2', 'Rsq_ML_s2', 'phase_ML_s2', 'amp_ML_s2', 'fitX_ML_s2', 'fitY_ML_s2',...
%         'A_HH_s0', 'B_HH_s0', 'phi_HH_s0', 'y_HH_s0', 'SSres_HH_s0', 'SStot_HH_s0', 'Rsq_HH_s0', 'phase_HH_s0', 'amp_HH_s0', 'fitX_HH_s0', 'fitY_HH_s0',...
%         'A_HL_s0', 'B_HL_s0', 'phi_HL_s0', 'y_HL_s0', 'SSres_HL_s0', 'SStot_HL_s0', 'Rsq_HL_s0', 'phase_HL_s0', 'amp_HL_s0', 'fitX_HL_s0', 'fitY_HL_s0',...
%         'A_MH_s0', 'B_MH_s0', 'phi_MH_s0', 'y_MH_s0', 'SSres_MH_s0', 'SStot_MH_s0', 'Rsq_MH_s0', 'phase_MH_s0', 'amp_MH_s0', 'fitX_MH_s0', 'fitY_MH_s0',...
%         'A_ML_s0', 'B_ML_s0', 'phi_ML_s0', 'y_ML_s0', 'SSres_ML_s0', 'SStot_ML_s0', 'Rsq_ML_s0', 'phase_ML_s0', 'amp_ML_s0', 'fitX_ML_s0', 'fitY_ML_s0');

    %% polar and linear plots with actual delta IBI amplitudes: linear plot also includes the sinusoidal fit!!!
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
    scatter(phase, amp, 20, 'k', 'filled')
    hold on
    plot(fitX,fitY, 'r-', 'LineWidth', 4)
    hold on
    plot(fitX,fitY-B,'r--','LineWidth',4)
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

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% s-3 ALL PLOT
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_s3, amp_s3);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI - s3): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    
    % Plot as polar scatter
    subplot(1,2,1)
    polarplot(phase_s3, amp_s3, 'ko')
    title('Polar Plot of Phase vs Delta IBI')
    set(gca, 'FontSize', 14)
    
    % Plot as linear scatter
    subplot(1,2,2)
    scatter(phase_s3, amp_s3, 40, 'k', 'filled')
    hold on
    plot(fitX_s3, fitY_s3, 'r-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    title([sprintf('Circular-linear Correlation\nrho = %.2f, p = %.4f', rho, pval) ' Subj: ' subjectNo]);
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Real Data-- Fitted Amplitude: %.3f | Phase Shift: %.1f° | Offset: %.3f | R-sq: %.3f', ...
        A_s3, rad2deg(phi_s3), B_s3, Rsq_s3);
    
    % Add it as the sgtitle
    sgtitle(titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s3_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s3_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% s-2 ALL PLOT
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_s2, amp_s2);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI - s2): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    
    % Plot as polar scatter
    subplot(1,2,1)
    polarplot(phase_s2, amp_s2, 'ko')
    title('Polar Plot of Phase vs Delta IBI')
    set(gca, 'FontSize', 14)
    
    % Plot as linear scatter
    subplot(1,2,2)
    scatter(phase_s2, amp_s2, 20, 'k', 'filled')
    hold on
    plot(fitX_s2, fitY_s2, 'r-', 'LineWidth', 4)
    hold on
    plot(fitX_s2, fitY_s2-B_s2, 'r--', 'LineWidth', 4)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    title([sprintf('Circular-linear Correlation\nrho = %.2f, p = %.4f', rho, pval) ' Subj: ' subjectNo]);
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Real Data-- Fitted Amplitude: %.3f | Phase Shift: %.1f° | Offset: %.3f | R-sq: %.3f', ...
        A_s2, rad2deg(phi_s2), B_s2, Rsq_s2);
    
    % Add it as the sgtitle
    sgtitle(titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s2_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s2_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% s-1 ALL PLOT
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_s1, amp_s1);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI - s1): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    
    % Plot as polar scatter
    subplot(1,2,1)
    polarplot(phase_s1, amp_s1, 'ko')
    title('Polar Plot of Phase vs Delta IBI')
    set(gca, 'FontSize', 14)
    
    % Plot as linear scatter
    subplot(1,2,2)
    scatter(phase_s1, amp_s1, 20, 'k', 'filled')
    hold on
    plot(fitX_s1, fitY_s1, 'r-', 'LineWidth', 4)
    hold on
    plot(fitX_s1, fitY_s1-B_s1, 'r--', 'LineWidth', 4)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    title([sprintf('Circular-linear Correlation\nrho = %.2f, p = %.4f', rho, pval) ' Subj: ' subjectNo]);
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Real Data-- Fitted Amplitude: %.3f | Phase Shift: %.1f° | Offset: %.3f | R-sq: %.3f', ...
        A_s1, rad2deg(phi_s1), B_s1, Rsq_s1);
    
    % Add it as the sgtitle
    sgtitle(titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s1_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s1_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% s-0 ALL PLOT
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_s0, amp_s0);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI - s0): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    
    % Plot as polar scatter
    subplot(1,2,1)
    polarplot(phase_s0, amp_s0, 'ko')
    title('Polar Plot of Phase vs Delta IBI')
    set(gca, 'FontSize', 14)
    
    % Plot as linear scatter
    subplot(1,2,2)
    scatter(phase_s0, amp_s0, 20, 'k', 'filled')
    hold on
    plot(fitX_s0, fitY_s0, 'r-', 'LineWidth', 4)
    hold on
    plot(fitX_s0, fitY_s0-B_s0, 'r--', 'LineWidth', 4)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    title([sprintf('Circular-linear Correlation\nrho = %.2f, p = %.4f', rho, pval) ' Subj: ' subjectNo]);
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Real Data-- Fitted Amplitude: %.3f | Phase Shift: %.1f° | Offset: %.3f | R-sq: %.3f', ...
        A_s0, rad2deg(phi_s0), B_s0, Rsq_s0);
    
    % Add it as the sgtitle
    sgtitle(titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s0_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_s0_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% s-p1 ALL PLOT
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_s_p1, amp_s_p1);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI - s-p1): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    
    % Plot as polar scatter
    subplot(1,2,1)
    polarplot(phase_s_p1, amp_s_p1, 'ko')
    title('Polar Plot of Phase vs Delta IBI')
    set(gca, 'FontSize', 14)
    
    % Plot as linear scatter
    subplot(1,2,2)
    scatter(phase_s_p1, amp_s_p1, 40, 'k', 'filled')
    hold on
    plot(fitX_s_p1, fitY_s_p1, 'r-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    title([sprintf('Circular-linear Correlation\nrho = %.2f, p = %.4f', rho, pval) ' Subj: ' subjectNo]);
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Real Data-- Fitted Amplitude: %.3f | Phase Shift: %.1f° | Offset: %.3f | R-sq: %.3f', ...
        A_s_p1, rad2deg(phi_s_p1), B_s_p1, Rsq_s_p1);
    
    % Add it as the sgtitle
    sgtitle(titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_sp1_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_sp1_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% linear plots with actual delta IBI amplitudes: Hit vs Miss!!!
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_hit, amp_hit);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI hit): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    scatter(phase_hit, amp_hit, 40, 'g', 'filled')
    hold on
    plot(fitX_hit,fitY_hit, 'g-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 14)

    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_miss, amp_miss);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI hit): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    hold on
    scatter(phase_miss, amp_miss, 40, 'r', 'filled')
    hold on
    plot(fitX_miss, fitY_miss, 'r-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Offset Hits: %.3f | R-sq Hits: %.3f | Offset Misses: %.3f | R-sq Misses: %.3f', ...
        B_hit, Rsq_hit, B_miss, Rsq_miss);
    title([titleStr ' | Subj: ' subjectNo]);
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_hitMiss_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_hitMiss_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% linear plots with actual delta IBI amplitudes: Fast vs Slow!!!
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_fast, amp_fast);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI fast): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    scatter(phase_fast, amp_fast, 40, 'g', 'filled')
    hold on
    plot(fitX_fast, fitY_fast, 'g-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 14)

    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_slow, amp_slow);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI slow): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    hold on
    scatter(phase_slow, amp_slow, 40, 'r', 'filled')
    hold on
    plot(fitX_slow, fitY_slow, 'r-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Offset Fast: %.3f | R-sq Fast: %.3f | Offset Slow: %.3f | R-sq Slow: %.3f', ...
        B_fast, Rsq_fast, B_slow, Rsq_slow);
    title([titleStr ' | Subj: ' subjectNo]);

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_fastSlow_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_fastSlow_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% linear plots with actual delta IBI amplitudes: highC vs lowC!!!
    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_highC, amp_highC);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI highC): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    figure;
    scatter(phase_highC, amp_highC, 40, 'g', 'filled')
    hold on
    plot(fitX_highC,fitY_highC, 'g-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 14)

    % Run circular-linear correlation
    [rho, pval] = circ_corrcl(phase_lowC, amp_lowC);
    
    % Display results
    fprintf('Circular-linear correlation (actual delta IBI highC): rho = %.3f, p = %.4f\n', rho, pval);
    
    % Visualization
    hold on
    scatter(phase_lowC, amp_lowC, 40, 'r', 'filled')
    hold on
    plot(fitX_lowC, fitY_lowC, 'r-', 'LineWidth', 2)
    xlabel('Breathing Phase (rad)')
    ylabel('Delta IBI')
    set(gca, 'FontSize', 14)

    % Format the string
    titleStr = sprintf('Offset highCs: %.3f | R-sq highCs: %.3f | Offset lowCs: %.3f | R-sq lowCs: %.3f', ...
        B_highC, Rsq_highC, B_lowC, Rsq_lowC);
    title([titleStr ' | Subj: ' subjectNo]);
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_highClowC_subj' subjectNo '.fig'];
    saveas(gcf,filename)
    filename = [saveDir '/circ2linear_corr_brPhase_ACD_highClowC_subj' subjectNo '.png'];
    saveas(gcf,filename)
    close;

    %% convert angles to radians
%     addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_statisticsScripts\circstat-matlab-master\circstat-matlab-master
%     min3_rad = circ_ang2rad(s_min3_allPhases);
%     min2_rad = circ_ang2rad(s_min2_allPhases);
%     min1_rad = circ_ang2rad(s_min1_allPhases);
%     zero_rad = circ_ang2rad(s_minZ_allPhases);
%     stim_rad = circ_ang2rad(s_min_stim_allPhases);
%     plus1_rad = circ_ang2rad(s_plus1_allPhases);
%     plus2_rad = circ_ang2rad(s_plus2_allPhases);
%     
%     min3_rad(isnan(min3_rad)) = [];
%     mu_min3 = circ_mean(min3_rad);
%     len_min3 = circ_r(min3_rad);
%     subplot(2,4,1)
%     ax_min3 = circ_plot(min3_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-3')
% 
%     min2_rad(isnan(min2_rad)) = [];
%     mu_min2 = circ_mean(min2_rad);
%     len_min2 = circ_r(min2_rad);
%     subplot(2,4,2)
%     ax_min2 = circ_plot(min2_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-2')
% 
%     min1_rad(isnan(min1_rad)) = [];
%     mu_min1 = circ_mean(min1_rad);
%     len_min1 = circ_r(min1_rad);
%     subplot(2,4,3)
%     ax_min1 = circ_plot(min1_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-1')
% 
%     zero_rad(isnan(zero_rad)) = [];
%     mu_zero = circ_mean(zero_rad);
%     len_zero = circ_r(zero_rad);
%     subplot(2,4,4)
%     ax_zero = circ_plot(zero_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-Zero')
% 
%     stim_rad(isnan(stim_rad)) = [];
%     mu_stim = circ_mean(stim_rad);
%     len_stim = circ_r(stim_rad);
%     subplot(2,4,5)
%     ax_stim = circ_plot(stim_rad,'pretty','ro',true,'linewidth',3,'color',[0.5 0 0]);
%     title('R-peak: Stim Onset')
% 
%     plus1_rad(isnan(plus1_rad)) = [];
%     mu_plus1 = circ_mean(plus1_rad);
%     len_plus1 = circ_r(plus1_rad);
%     subplot(2,4,6)
%     ax_plus1 = circ_plot(plus1_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-Plus 1')
%     
%     plus2_rad(isnan(plus2_rad)) = [];
%     mu_plus2 = circ_mean(plus2_rad);
%     len_plus2 = circ_r(plus2_rad);
%     subplot(2,4,7)
%     ax_plus2 = circ_plot(plus2_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-Plus 2')
% 
%     screenSize = get(0, 'ScreenSize');
%     set(gcf, 'Position', screenSize);
% 
%     saveDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo];
%     figFileName = 'allTrials_breathingPhaseLocking.fig';
%     % Save the figure
%     saveas(gcf, fullfile(saveDir, figFileName));
%     figFileName = 'allTrials_breathingPhaseLocking.png';
%     saveas(gcf, fullfile(saveDir, figFileName));
%     close;
% 
%     fileName = 'breathing_vectorLengths_meanAngles.mat';
% 
%     save(fullfile(saveDir, fileName), 'mu_min3', 'len_min3', 'mu_min2', 'len_min2', 'mu_min1', 'len_min1', 'mu_zero', 'len_zero', 'mu_stim', 'len_stim','mu_plus1','len_plus1','mu_plus2','len_plus2');
    
    %% convert angles to radians in the exhale and inhale trial phases separately
%     addpath S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_statisticsScripts\circstat-matlab-master\circstat-matlab-master
% %     min3_rad = circ_ang2rad(s_min3_allPhases);
%     min2_exh_rad = circ_ang2rad(exhale_s2Phase);
%     min1_exh_rad = circ_ang2rad(exhale_s1Phase);
%     zero_exh_rad = circ_ang2rad(exhale_s0Phase);
%     min2_inh_rad = circ_ang2rad(inhale_s2Phase);
%     min1_inh_rad = circ_ang2rad(inhale_s1Phase);
%     zero_inh_rad = circ_ang2rad(inhale_s0Phase);
%     
%     min2_exh_rad(isnan(min2_exh_rad)) = [];
%     mu_min2_exh = circ_mean(min2_exh_rad);
%     len_min2_exh = circ_r(min2_exh_rad);
%     subplot(2,3,1)
%     ax_min2_exh = circ_plot(min2_exh_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.3]);
%     title('R-peak: Stim-2 Exhale')
%     
%     min2_inh_rad(isnan(min2_inh_rad)) = [];
%     mu_min2_inh = circ_mean(min2_inh_rad);
%     len_min2_inh = circ_r(min2_inh_rad);
%     subplot(2,3,4)
%     ax_min2_inh = circ_plot(min2_inh_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.3]);
%     title('R-peak: Stim-2 Inhale')
% 
%     min1_exh_rad(isnan(min1_exh_rad)) = [];
%     mu_min1_exh = circ_mean(min1_exh_rad);
%     len_min1_exh = circ_r(min1_exh_rad);
%     subplot(2,3,2)
%     ax_min1_exh = circ_plot(min1_exh_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-1 Exhale')
%     
%     min1_inh_rad(isnan(min1_inh_rad)) = [];
%     mu_min1_inh = circ_mean(min1_inh_rad);
%     len_min1_inh = circ_r(min1_inh_rad);
%     subplot(2,3,5)
%     ax_min1_inh = circ_plot(min1_inh_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.5]);
%     title('R-peak: Stim-1 Inhale')
% 
%     zero_exh_rad(isnan(zero_exh_rad)) = [];
%     mu_zero_exh = circ_mean(zero_exh_rad);
%     len_zero_exh = circ_r(zero_exh_rad);
%     subplot(2,3,3)
%     ax_zero_exh = circ_plot(zero_exh_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.7]);
%     title('R-peak: Zero Exhale')
%     
%     zero_inh_rad(isnan(zero_inh_rad)) = [];
%     mu_zero_inh = circ_mean(zero_inh_rad);
%     len_zero_inh = circ_r(zero_inh_rad);
%     subplot(2,3,6)
%     ax_zero_inh = circ_plot(zero_inh_rad,'pretty','bo',true,'linewidth',3,'color',[0 0 0.7]);
%     title('R-peak: Zero Inhale')
% 
%     screenSize = get(0, 'ScreenSize');
%     set(gcf, 'Position', screenSize);
% 
%     saveDir = ['S:\KNEU\KNEUR-Projects\Projects\Ege-Backup\VTD\sukanya_MScThesis\Sukanya-Backup\VTD_CircularStats\SubjectwiseData\VTD' subjectNo];
%     figFileName = 'breathingPhaseLocking_inhaleVsExhale.fig';
%     % Save the figure
%     saveas(gcf, fullfile(saveDir, figFileName));
%     figFileName = 'breathingPhaseLocking_inhaleVsExhale.png';
%     saveas(gcf, fullfile(saveDir, figFileName));
%     close;
% 
%     fileName = 'breathing_vectorLengths_meanAngles_inhaleVsExhale.mat';
% 
%     save(fullfile(saveDir, fileName), 'mu_min2_exh', 'len_min2_exh', 'mu_min2_inh', 'len_min2_inh', 'mu_min1_exh', 'len_min1_exh', 'mu_min1_inh', 'len_min1_inh','mu_zero_exh', 'len_zero_exh', 'mu_zero_inh', 'len_zero_inh');
    
%% categorize trials according to fast and slow-response trials
% medianRespTime = median(allRespTimes,'omitnan');
% fastResps = [];
% slowResps = [];
% if size(allRespTimes,2) == size(breath_set.epoch,2)
%     for tr=1:size(breath_set.epoch,2)
%         if allRespTimes(tr) <= medianRespTime
%             fastResps = [fastResps tr];
%         elseif allRespTimes(tr) > medianRespTime
%             slowResps = [slowResps tr];
%         end
%     end
% else
%     disp('ERROR in trial sorting, check Epoch Events!!!')
%     return
% end
% filename_resp = 'fast_slowResponseTrials.mat';
%     save(fullfile(saveDir, filename_resp), 'fastResps','slowResps');

%% categorize trials according to low and high confidence
%     array1 = confTable.trialIBIs(:,1);
%     array2 = ibi_s1';
%     
%     n = length(array1);
%     m = length(array2);
%     
%     % Initialize LCS table
%     L = zeros(n+1, m+1);
%     
%     % Build the table
%     for i = 1:n
%         for j = 1:m
%             if array1(i) == array2(j)
%                 L(i+1, j+1) = L(i, j) + 1;
%             else
%                 L(i+1, j+1) = max(L(i+1, j), L(i, j+1));
%             end
%         end
%     end
%     
%     % Backtrack to find indices of matching elements
%     i = n; j = m;
%     matchIdx1 = [];
%     matchIdx2 = [];
%     
%     while i > 0 && j > 0
%         if array1(i) == array2(j)
%             matchIdx1(end+1) = i;
%             matchIdx2(end+1) = j; 
%             i = i - 1;
%             j = j - 1;
%         elseif L(i, j+1) > L(i+1, j)
%             i = i - 1;
%         else
%             j = j - 1;
%         end
%     end
%     
%     matchIdx1 = sort(matchIdx1);
%     matchIdx2 = sort(matchIdx2);
%     
%     % Find extra elements
%     allIdx1 = (1:n)';
%     allIdx2 = (1:m)';
%     extraIdx1 = setdiff(allIdx1, matchIdx1);
%     extraIdx2 = setdiff(allIdx2, matchIdx2);
%     
%     if isempty(extraIdx2)
%         % Optionally remove them
%         array1_cleaned = array1;
%         array1_cleaned(extraIdx1) = [];
%         % Display result
%         disp('Indices of extra elements in array1:');
%         disp(extraIdx1);
%         
%         confClean = confTable(setdiff(1:size(confTable,1),extraIdx1),:);
%         highConfTrials = find(confClean.confidence>medianConf);
%         lowConfTrials = find(confClean.confidence<=medianConf);
%         filename_conf = 'lowHighConfTrials.mat';
%         save(fullfile(saveDir, filename_conf), 'highConfTrials','lowConfTrials');
%     else
%         array2(extraIdx2) = [];
%         disp('some cleaning was necessary');
%         return
%     end



end
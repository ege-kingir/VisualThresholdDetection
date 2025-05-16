function [output, keyboard_events, params, ET, gazes] = trials_VTDnew_egeVersion(params)
% Trials function for MAIN_VTD.m

% Sukanya C, Ege Kingir (MSc Thesis, IMPRS)
% 12 Nov 2023
% Trials function for the visual detection task after threshold calibration

%%% This function contains the trial loop. The loop is a state-based one:
%%% If the subject fixates to the center successfully (state=1), we move on
%%% to the warning cue state (state=2). Afterwards, the stimulus is
%%% presented. Then some delay and "detection (yes/no)" and "confidence"
%%% reports by the subject.

%%% If the subject breaks fixation to the center during the pre-stimulus
%%% window, the trial restarts.

%%% 60 trials per block: 87% stimulus-present, 13% catch.
%%% 6-8 blocks in a session.

%% CHANGE LOG
% Following changes made by XX on XX.XX.XXX
% 30.11.2023 -- EK updates:
% 1) Added the continuous confidence scale bar.

% 2) Under state==1, flag2 should be defined as zero before you enter the...
% %...loop. Otherwise the variable is not recognized (I added the...
% %...definition).

% 13.12.2023 updates:
% Removed the online plots for gaze and accuracy as they are slowing the
% code down 


%% FUNCTIONS CALLED IN THIS SCRIPT / RELEVANT FILES
% params_VTD.m
% startupEyetracking.m
% SendSignal.m
% display_cycle_clock.m (only if the task involves cycling)


%% Set up eye tracking
if params.eyetracking == 1
    if params.whicheyetracker == 1 %Eyelink eyetracker (currently not in use)
        params.samplingrate = 1000;

        params.maxeyesamples = params.numtrials*params.totaltrialdur*params.samplingrate;

        params.eyedatafile = sprintf('M%s_%d.edf', params.subject, params.runnumber); % has to be less than 8 characters

        % Initialize Eyelink
        el = startupEyetracking(params.window, params.eyedatafile);

        % start recording eye position
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        WaitSecs(0.1);
        % mark zero-plot time in data file
        Eyelink('Message', 'SYNCTIME');

        % Initialize eye data structure
        eyetrack.trial(1:params.maxeyesamples,1) = NaN;
        eyetrack.state(1:params.maxeyesamples,1) = NaN;
        eyetrack.eyetime(1:params.maxeyesamples,1) = NaN;
        eyetrack.ptbtime(1:params.maxeyesamples,1) = NaN;
        eyetrack.xpos(1:params.maxeyesamples,2) = NaN;
        eyetrack.ypos(1:params.maxeyesamples,2) = NaN;
        eyetrack.pupil(1:params.maxeyesamples,2) = NaN;
        
    elseif params.whicheyetracker == 2

         %%% Initialize Pupil Labs Remote Control
            % Pupil Remote address -- you can check this from PupilCapture
            endpoint =  'tcp://169.254.53.40:50020';

            % Setup zmq context and remote helper
            ctx = zmq.core.ctx_new();
            socket = zmq.core.socket(ctx, 'ZMQ_REQ');
            % set timeout to 1000ms in order to not get stuck in a blocking
            % mex-call if server is not reachable, see
            % http://api.zeromq.org/4-0:zmq-setsockopt#toc19
            zmq.core.setsockopt(socket, 'ZMQ_RCVTIMEO', 1000);

            fprintf('Connecting to %s\n', endpoint);
            zmq.core.connect(socket, endpoint);
            
            % Request sub port
            zmq.core.send(socket, uint8('SUB_PORT'));
            sub_port = char(zmq.core.recv(socket));
            fprintf('Received sub port: %s\n', sub_port);
            
            if isequal(sub_port, false)
                warning('No valid sub port received');
                return;  % exit script
            end
            
            % Create and connect sub socket
            ip_address = '169.254.53.40';
            sub_endpoint =  sprintf('tcp://%s:%s', ip_address, sub_port);
            sub_socket = zmq.core.socket(ctx, 'ZMQ_SUB');
            %Connect to sub port
            zmq.core.connect(sub_socket, sub_endpoint);
            zmq.core.setsockopt(sub_socket, 'ZMQ_RCVTIMEO', 1000); 
            fprintf('Connecting to SUB: %s\n', sub_endpoint);
                        
            tic; % Measure round trip delay
            zmq.core.send(socket, uint8('t'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));
            fprintf('Round trip command delay: %s\n', toc);

            % set current Pupil time to 0.0
            zmq.core.send(socket, uint8('T 0.0'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));

            % start recording
            pause(1.0);
            zmq.core.send(socket, uint8('R')); %%%%
            result = zmq.core.recv(socket); %%%%
            fprintf('Recording should start: %s\n', char(result));

        %     zmq.core.send(socket, uint8('t'));
        %     currentTime = zmq.core.recv(socket);
    end;
    
end;

%% Pupil Labs calibration
if params.eyetracking == 1
    
    if params.whicheyetracker == 2
        
        %open a grey window
        [params.window, params.windowRect] = PsychImaging('OpenWindow', params.screenNumber, params.grey);
        
        % Size of the on screen window in pixels
        %     [params.screenXpixelscalib, params.screenYpixelscalib] = Screen('WindowSize', params.window);
        Screen('BlendFunction', params.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
        
        %get the center coordinates (in pixels, bottom-left of the screen
        %is [0 0].
        [params.screenXcentre_pix, params.screenYcentre_pix] = RectCenter(params.windowRect);
        params.centercoords = [params.screenXcentre_pix params.screenYcentre_pix];
        %% get the pupil labs calibration marker
        % % marker = imread('calibMarker.jpg');
        % % markerTexture = Screen('MakeTexture', params.window, marker);
        % %
        % % Screen('DrawTexture', params.window, markerTexture, [], [], 0);
        baseRect1 = [0 0 165 165];
        baseRect2 = [0 0 160 160];
        baseRect3 = [0 0 100 100];
        baseRect4 = [0 0 60 60];
        baseRect5 = [0 0 30 30];
        x=[0 -params.screenX_pix/2+baseRect1(3) params.screenX_pix/2-baseRect1(3) -params.screenX_pix/2+baseRect1(3) params.screenX_pix/2-baseRect1(3)];
        y=[0 -params.screenY_pix/2+baseRect1(3) -params.screenY_pix/2+baseRect1(3) params.screenY_pix/2-baseRect1(3) params.screenY_pix/2-baseRect1(3)];
        valid = 1;
        for valid=1:4
            for stim=1:5
                
                centeredRect1 = CenterRectOnPointd(baseRect1, params.screenXcentre_pix+x(stim), params.screenYcentre_pix+y(stim));
                centeredRect2 = CenterRectOnPointd(baseRect2, params.screenXcentre_pix+x(stim), params.screenYcentre_pix+y(stim));
                centeredRect3 = CenterRectOnPointd(baseRect3, params.screenXcentre_pix+x(stim), params.screenYcentre_pix+y(stim));
                centeredRect4 = CenterRectOnPointd(baseRect4, params.screenXcentre_pix+x(stim), params.screenYcentre_pix+y(stim));
                centeredRect5 = CenterRectOnPointd(baseRect5, params.screenXcentre_pix+x(stim), params.screenYcentre_pix+y(stim));
                
                xCoords = [-2 2 0 0];
                yCoords = [0 0 -2 2];
                allCoords = [xCoords; yCoords];
                lineWidthPix = 2;
                
                %nested doughnut shaped ovals to make calibration spots...
                %recognizable by PupilCapture.
                Screen('FillOval', params.window, [0.5 0.5 0.5],centeredRect1, max(baseRect1)*1.01);
                Screen('FillOval', params.window, [1 1 1],centeredRect2, max(baseRect2)*1.01);
                Screen('FillOval', params.window, [0 0 0],centeredRect3, max(baseRect3)*1.01);
                Screen('FillOval', params.window, [1 1 1],centeredRect4, max(baseRect4)*1.01);
                Screen('FillOval', params.window, [0 0 0],centeredRect5, max(baseRect5)*1.01);
                
                %fixation cross
                Screen('DrawLines', params.window, allCoords, lineWidthPix, [1 1 1], [params.screenXcentre_pix+x(stim) params.screenYcentre_pix+y(stim)],2);
                % Screen('DrawLines', params.window,
                Screen('Flip', params.window);
                KbStrokeWait;
                
            end
            Screen('Flip', params.window);
%             calibstart = input('Do you want to restart the calibration (in case accuracy is low) (1 = Yes, 0 = No) ');
%             if calibstart == 1
%                 valid = 1;
%             elseif calibstart == 0
%                 valid = 0;
%                 %
%             end
            
        end
        WaitSecs(5);
        Screen('Flip',params.window);
    end
end

%% Initialisation of parameters for checking fixation
ET.TrialAborted = false;
ET.FixationEstablished = false;
ET.FixationPending = false;
ET.FixationPendingStartTime = NaN;
ET.AbortPending = false;
ET.AbortPendingStartTime = NaN;
ET.AbortStartTime = NaN;
ET.IsFixation = false;

%% Cycling warm up if necessary
if strcmp(params.runtype, 'TC') || strcmp(params.runtype, 'L') || strcmp(params.runtype, 'H')
    display_cycle_clock_VTD(params);
end


%% Output struct
output = struct('condition', [], 'trialonset', []);


%% Initialising 

triggertime = NaN;
triggerout = 0;
cursample = 0; 


%% Keyboard set up
% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems:
KbName('UnifyKeyNames');

% Prevent spilling of keystrokes into console:
ListenChar(-1); % suppress keyboard input

% Define key list
keylist = zeros(1,256); % create mask for keys
% set keylist / mask to 1 for all the keys you want to monitor
% keylist(KbName('1!')) = 1; 
keylist(KbName('3#')) = 1;
keylist(KbName('2@')) = 1;
keylist(KbName('1!')) = 1;
% keylist(KbName('5%')) = 1;

% create and start a restricted KbQueue
KbQueueCreate(-1, keylist); %'-1' uses default keyboard, this is a restricted KbQueue, used for detecting button presses
KbQueueStart;

% Initialize keyboard_events as an empty structure array
keyboard_events = struct([]);

trial = 0;
flag = 0;

% Define your response keys
responseKeys = {'3#', '2@', '1!'};


%% Set up fixed and variable delays

fixation_Duration = 2.7; %% In Park et al., fixation duration is specified as 0.5 to 0.7 seconds, 
% and if stable gaze fixation 
% is achieved for the said duration, a trial is initiated. However, since
% it is yet not well sorted how to feed in Pupil Capture data into MATLAB,
% we choose to define a range of 0.7s (max gaze fixation defined in Park et
% al., 2014) + 0.5 seconds i.e., 1.2 seconds. 

warning_before_stimRange = [1.5, 1.7]; % this is the variable delay for which the fixation remains red before the stimulus appears
% stimulusDuration = 0.05; % this is the fixed stimulus presentation duration in seconds, converted to frames
stimulusFrames = round(0.05/params.flipinterval);
stimulusDuration = 0.05;
var_delay_after_stimRange = [0.3, 0.7]; % this is the variable delay range after stimulus presentation
time_for_response = 3; % this is the time accounted to the participant to make a key press

%% Initialise parameters for the QUEST algorithm 

tGuess = 0.04; % Estimated signal strength (contrast value) which is the threshold for detection
tGuessSd = 2; % Standard Deviation allowed to the QUEST algorithm for the distribution centred on tGuess
pThreshold = 0.65; % Probability of correct response - We want the hit rate to be at 65%
beta = 3.5; % Slope of the Weibull function
delta = 0.01;
gamma = 0; % Guess rate or probabilty of an incorrect guess at values high above threshold
grain = 0.0025; % Step shift in contrast - WE DO NOT CHANGE THE GRAIN
range = 0.075; %Range of contrast values under practical testing conditions
keyHit = '2@'; % Key for hit
keyMiss = '3#'; % Key for miss

% Initialize QUEST parameters
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma, grain, range);
p = QuestPdf(q,tGuess);
% q.x = linspace(0,1,501);
% q.x2 = linspace(0,2,1001);

% Define the center of the screen as empty!
center = [];

% Define the iterations fo getting gaze on surface inputs
iter=1;

for block = 1:params.numblocks
    
    random_trl_order = params.trialorder;
    Screen('TextSize',params.window, 30)
    Screen('DrawText',params.window, 'FIXATION CALIBRATION!',(params.screenX_pix/2)-100, params.screenY_pix*0.5, params.black);
    Screen('Flip', params.window);
    WaitSecs(1);
 
    % Plotting parameters initialisation
    hitCount = 0;
    missCount = 0;
    faCount = 0;
    hitRates = zeros(1, length(random_trl_order)); % Initialize hit rates vector
    falseAlarmRates = zeros(1, length(random_trl_order)); % Initialize false alarm rates vector
    stimPresentCount = 0;
    stimAbsentCount = 0;
    contrasts = zeros(1, length(random_trl_order)); % Initialize contrasts vector

    x1 = [];
    y1 = [];
    x2 = [];
    y2 = [];
    x3 = [];
    y3 = [];
    
%     figure(1);
%     subplot(2, 1, 1);
%     % % % plot1 = plot(x1, y1, 'r-');  % Initialize plot 1 (red line)
%     xlabel('Trial');
%     ylabel('Hit and False Alarm Rates');
%     title('QUEST modulated accuracy');

    
%     subplot(2, 1, 2);
%     % % % plot3 = plot(x3, y3, 'k-');  % Initialize plot 1 (red line)
%     xlabel('Trial');
%     ylabel('Grating Contrast');
%     title('QUEST modulated contrast');

   
    %figure(10);
    % % % plot1 = plot(x1, y1, 'r-');  % Initialize plot 1 (red line)
%     xlabel('x gaze');
%     ylabel('y gaze');
%     title('Real Time gaze');


%% Main Trials Loop
    %% Display April tags on PTB screen to detect surface via PupilLabs
    Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
    Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
    Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
    Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
    % fixation circle
    Screen('FillOval', params.window, params.black, params.centeredcircleRectOUT, params.maxDiameterOUT);
    Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
    Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, params.white, [], 1);
    Screen('Flip', params.window);
    if params.eyetracking == 1
        if params.whicheyetracker == 2
            zmq.core.setsockopt(sub_socket, 'ZMQ_SUBSCRIBE', 'surface'); %now you can subscribe to surface because we drew it!
        end
    end
    WaitSecs(2);
    
    initialSurfaceDef = 3; %number of secs allowed to define the surface and at the same time get baseline gaze data
    startSurfDef = GetSecs();
    gazes = [];
    while 1
        Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
        Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
        Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
        Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
        % fixation circle
        Screen('FillOval', params.window, params.black, params.centeredcircleRectOUT, params.maxDiameterOUT);
        Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
        Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, [0.8 0.8 0.8], [], 1);
        Screen('Flip', params.window);
        if params.eyetracking == 1
            if params.whicheyetracker == 2
                cursample = cursample + 1;
                % Get gaze positions from pupil labs
                % (filtermessages)
                [topic, note] = recv_message(sub_socket, 2500);
                if ~isequal(note, false)  % test for valid message
                    gaze_positions{1,iter}= [topic, note('gaze_on_surfaces')];   % print pupil norm_pos
                    for j=2:size(gaze_positions{1,iter},2)
                        gazecoordinates = gaze_positions{1,iter}{1,j}('norm_pos');
                        x(j) = gazecoordinates{1,1};
                        y(j) = gazecoordinates{1,2};

                    end

                    %This alternative is decreasing the temporal...
                    %...resolution in case matlab fails.
                    ET.CurrentGaze(1,1) = x(j);
                    ET.CurrentGaze(1,2) = y(j);
                    %                         plot(x(j),y(j),'p')
                    %                         drawnow
                    %                         hold on
                end
                gazes(end+1,1) = ET.CurrentGaze(1,1);
                gazes(end,2) = ET.CurrentGaze(1,2);
                ET.currentTime = GetSecs();
                ET.eyetime(cursample,1) = ET.currentTime;
    %                 ET = CheckFixation_VTD_egeVersion(ET, trial, params, output);
                iter=iter+1;
            end

        end

        if GetSecs()> startSurfDef + initialSurfaceDef
            break
        end
    end
    meanGazeOne = [mean(gazes(:,1)) mean(gazes(:,2))];
    meanGazeOne
    %% start the trial loop here
    while 1
        
        trial=trial+1;
        if random_trl_order(trial) == 1 % For stimulus present trials
            stimPresentCount = stimPresentCount + 1; % Increment stimulus present trial count
            % WaitSecs(fixation_Duration);
        elseif random_trl_order(trial) == 2 % For stimulus absent trials
            stimAbsentCount = stimAbsentCount + 1; % Increment stimulus absent count
            % WaitSecs(fixation_Duration);
        end
        % We only have a single trial type for now, as we are not considering
        % the inter-trial interval separately, but rather adding it in with the
        % fixation check to initiate trial
        if flag == 1
            Screen('DrawText', params.window, 'FIXATE PLEASE!', (params.screenX_pix/2)-100, (params.screenY_pix/2), params.black);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            Screen('Flip', params.window);
            WaitSecs(2);
        end
        
        ET.TrialAborted = false;
        ET.FixationEstablished = false;
        ET.FixationPending = false;
        ET.FixationPendingStartTime = NaN;
        ET.AbortPending = false;
        ET.AbortPendingStartTime = NaN;
        ET.AbortStartTime = NaN;
        ET.IsFixation = false;
        flag = 0;
        earlyResponseFlag=0;
        flagResponse=1;

        %% Display durations that are jittered -- thus values need to change every trial
        output(trial).condition = random_trl_order(trial);
        delayBeforeStimulus = warning_before_stimRange(1) + rand() * (warning_before_stimRange(2) - warning_before_stimRange(1));
        delayAfterStimulus = var_delay_after_stimRange(1) + rand() * (var_delay_after_stimRange(2) - var_delay_after_stimRange(1));
        intertrial_Interval_Range = [1.5, 2]; % this is the variable delay range for intertrial interval in seconds - default in Park et al. ,
                                              % is 1.5 to 2, but since we allow another additional 1 sec to the subjects to start the fixation, 
                                              % we have reduced this range accordingly
        intertrialInterval = intertrial_Interval_Range(1) + rand() * (intertrial_Interval_Range(2) - intertrial_Interval_Range(1));
        params.TimeToBeginFixation = intertrialInterval; % Time allowed after the ITI to start fixating
        establish_fixation_range = [0.5, 0.7];
        establishFixation = establish_fixation_range(1) + rand() * (establish_fixation_range(2) - establish_fixation_range(1));
        params.TimeToEstablishFixation = establishFixation; % Time allowed to establish fixation after the time to begin fixation
        params.BlinkAllowance = 0.2; % Time given to allow for blinks, usually around 200-300 ms.
        
        
        %% target stimulus contrast in the current trial
        currentEstimate = QuestMean(q);
        output(trial).AllContrastEstimatesMean = currentEstimate;
        currentContrast = QuestQuantile(q);
        output(trial).AllContrastEstimatesQuantile = currentContrast;
        contrasts(trial) = currentContrast;
        
        %% STATE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nTrial %d of %d, Run %d, Trial Type - SP (1) or SA (2): %d \n\n\n\n\n\n', trial, params.numtrials, params.runnumber, random_trl_order(trial));
                
        % Hereon starts the first state - fixation on the central black circle
        state = 1; 
        flag2 = 0;
                   
        Screen('FillOval', params.window, params.black, params.centeredcircleRectOUT, params.maxDiameterOUT);
        Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
        Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, params.black, [], 1);
        %% Display April tags on PTB screen to detect surface via PupilLabs
        Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
        Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
        Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
        Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
        vbl = Screen('Flip', params.window);
        output(trial).trialonset = vbl;
        
        % Send EEG trigger code for trial start
        if params.EEG == 1
            SendSignal(params.ioObj, params.address, random_trl_order(trial)); %1 for stimulus present, 2 for stimulus absent
            triggerout = 1;
            triggertime = GetSecs;
        else
        end
        
        if params.eyetracking == 1 % Send Annotation from Matlab to Pupil Labs
            if params.whicheyetracker == 2
                zmq.core.send(socket, uint8('t'));
                currentTime_bytes = zmq.core.recv(socket,20);
                % fprintf('first measure: %s\n', char(currentTime_bytes));
                currentTimeChar = char(currentTime_bytes);
                currentTimeNum = str2num(currentTimeChar);
                keys_start = {'topic', 'label', 'timestamp', 'duration'};
                values_start = {'annotation.TrialStart', strcat('Trial #', num2str(trial),' started'), currentTimeNum, 1.0};
                start_annotation = containers.Map(keys_start, values_start);
                send_annotation(socket, start_annotation);
                result = zmq.core.recv(socket);
            end
        end
        
        gazes = []; %to note down all of the recorded gaze locations within the surface, during a trial! It is redefined as empty when the trial restarts due to fixation.
        if params.checkfix == 1
            if trial == 1
                center(trial,:) = meanGazeOne
                %                                 center
            elseif trial == 2
                oldGaze = meanGazeOne * params.weightOld;
                newGaze = output(1).meanGazes * params.weightNew;
%                 vector = [meanGazeOne;output(1).meanGazes];
                center(trial,:) = oldGaze + newGaze
                %                                 center
            else
                oldGaze = center(trial-1,:) * params.weightOld;
                newGaze = output(trial-1).meanGazes * params.weightNew;
%                 vector = [center(trial-1,:);output(trial-1).meanGazes];
                center(trial,:) = oldGaze + newGaze
                %                                 center
            end
        end
        
        while state == 1 

            [pressed, firstPress] = KbQueueCheck;
            
            if pressed
                % Extract information about the key presses
                [key,ind] = max(firstPress);
                % Check if any of the response keys are pressed

                if ismember(ind, KbName(responseKeys(1:2)))
                    % Identify the pressed key
                    % idx = find(ismember(KbName(responseKeys), pressedKeys));
                    % response = responseKeys{idx}; % Assign the pressed key as response
                    Screen('TextSize',params.window, 30);
                    Screen('DrawText',params.window, 'WAIT QUESTION MARK TO RESPOND!',(params.screenX_pix/2)-100, params.screenY_pix*0.5, params.black);
%                     Screen('DrawText',params.window, 'Please wait till the question mark sign appears on screen to make your report.',(params.screenX_pix/2)-100, params.screenY_pix*0.55, params.black);
                    %% Display April tags on PTB screen to detect surface via PupilLabs
                    Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                    Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                    Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                    Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                    Screen('Flip', params.window);
                    WaitSecs(2);
                    flag2 = 1;
                    
                end
               
            end
            
            if flag2 == 1
                fprintf('early response');
                earlyResponseFlag=1; %apart from trial abortion by fixation break; this is another way to abort the trial: because the subject pressed a response button earlier than they were supposed to!
                break
                
            end
            
            % If eyetracking
            %
            if params.eyetracking == 1
                if params.whicheyetracker == 2
                    if params.checkfix == 1
                        cursample = cursample + 1;
                        % Get gaze positions from pupil labs
                        % (filtermessages)
                        [topic, note] = recv_message(sub_socket, 2500);
                        if ~isequal(note, false)  % test for valid message
                            gaze_positions{1,iter}= [topic, note('gaze_on_surfaces')];  % print pupil norm_pos
                            for j=2:size(gaze_positions{1,iter},2)
                                gazecoordinates = gaze_positions{1,iter}{1,j}('norm_pos');
                                x(j) = gazecoordinates{1,1};
                                y(j) = gazecoordinates{1,2};
                                
                            end
                           
                        end
                        
                        ET.CurrentGaze(1,1) = x(j);
                        ET.CurrentGaze(1,2) = y(j);
                        
                        if GetSecs() > output(trial).trialonset + params.TimeToBeginFixation % we do not consider the gazes during the ITI for calculating the center for the next trial
                            gazes(end+1,1) = ET.CurrentGaze(1,1);
                            gazes(end,2) = ET.CurrentGaze(1,2);
                        end
                        
                        ET.currentTime = GetSecs();
                        ET.eyetime(cursample,1) = ET.currentTime;
                            %                                 
                        dx = abs(x(j) - center(trial,1));
                        dy = abs(y(j) - center(trial,2));

                        ET.offset_x = dx * params.pix_x;
                        ET.offset_y  = dy * params.pix_y;
                        ET = CheckFixation_VTD_egeVersion(ET, trial, params, output);

                    end
                end

            end
            
            if ET.TrialAborted == true
                flag = 1;
                if params.EEG == 1 
                     SendSignal(params.ioObj, params.address,55); %trial aborted trigger
                     triggerout = 1;
                     triggertime = GetSecs;
                else
                end
                %% Display April tags on PTB screen to detect surface via PupilLabs
                Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                Screen('DrawText', params.window, 'FIXATE PLEASE!', (params.screenX_pix/2)-100, (params.screenY_pix/2), params.black);
                Screen('Flip',params.window);
                break
            end
            
            if params.checkfix==1
                if ET.FixationEstablished && GetSecs>output(trial).trialonset+params.TimeToBeginFixation
                    state = 2;
                end
            else
                if GetSecs>output(trial).trialonset+params.TimeToBeginFixation
                    state = 2;
                end
            end
        end
        
        if flag == 1 || earlyResponseFlag ==1 %if trial is aborted, 'continue' comment takes you back to line with 'while 1'
            if random_trl_order(trial) == 1
                stimPresentCount=stimPresentCount-1;
            elseif random_trl_order(trial) == 2
                stimAbsentCount=stimAbsentCount-1;
            end
            trial = trial - 1;
            continue
        end
        
        %% STATE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if random_trl_order(trial) == 1 % For stimulus present trials
            
            % Change fixation circle to red and wait for
            % variable delay before stimulus presentation
            Screen('FillOval', params.window, [0.5 0 0], params.centeredcircleRectOUT, params.maxDiameterOUT);
            Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
            Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, [0.5 0 0], [], 1);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            Screen('Flip', params.window);
            % Get warning onset time
            output(trial).warningOnset = GetSecs;

        elseif random_trl_order(trial) == 2 % For stimulus absent trials
            
            % Change fixation circle to red and wait for
            % variable delay before stimulus presentation
            Screen('FillOval', params.window, [0.5 0 0], params.centeredcircleRectOUT, params.maxDiameterOUT);
            Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
            Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, [0.5 0 0], [], 1);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            Screen('Flip', params.window);
            % Get warning onset time
            output(trial).warningOnset = GetSecs;
        end
        
        flag3=0;
        earlyResponseFlag=0;
        while state == 2
            [pressed, firstPress] = KbQueueCheck;
            
            if pressed
                % Extract information about the key presses
                [~,pressedKeys] = find(firstPress);
                
                % Check if any of the response keys are pressed
                if any(ismember(pressedKeys, KbName(responseKeys(1:2))))
                    % Identify the pressed key
                    % idx = find(ismember(KbName(responseKeys), pressedKeys));
                    % response = responseKeys{idx}; % Assign the pressed key as response
                    Screen('TextSize',params.window, 30)
                    Screen('DrawText',params.window, 'WAIT QUESTION MARK TO RESPOND!',(params.screenX_pix/2)-100, params.screenY_pix*0.5, params.black);
%                     Screen('DrawText',params.window, 'Please wait till the question mark sign appears on screen to make your report.',(params.screenX_pix/2)-100, params.screenY_pix*0.55, params.black);
                    %% Display April tags on PTB screen to detect surface via PupilLabs
                    Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                    Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                    Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                    Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                    Screen('Flip', params.window);
                    WaitSecs(2);
                    flag3 = 1;
                    
                end
            end
            if flag3 == 1
                fprintf('early response');
                earlyResponseFlag = 1;
                break
            end
            % If eyetracking
            if params.eyetracking == 1
                if params.whicheyetracker == 2
                    if params.checkfix == 1
                        cursample = cursample + 1;
                        % Get gaze positions from pupil labs
                        % (filtermessages)
                        [topic, note] = recv_message(sub_socket, 2500);
                        if ~isequal(note, false)  % test for valid message
                            gaze_positions{1,iter}= [topic, note('gaze_on_surfaces')];  % print pupil norm_pos
                            for j=2:size(gaze_positions{1,iter},2)
                                gazecoordinates = gaze_positions{1,iter}{1,j}('norm_pos');
                                x(j) = gazecoordinates{1,1};
                                y(j) = gazecoordinates{1,2};
                                
                            end
                            %This alternative is decreasing the temporal...
                            %...resolution in case matlab fails.
                        end
                        
                        ET.CurrentGaze(1,1) = x(j);
                        ET.CurrentGaze(1,2) = y(j);
                        gazes(end+1,1) = ET.CurrentGaze(1,1);
                        gazes(end,2) = ET.CurrentGaze(1,2);
                        ET.currentTime = GetSecs();
                        ET.eyetime(cursample,1) = ET.currentTime;
                        dx = abs(x(j) - center(trial,1));
                        dy = abs(y(j) - center(trial,2));

                        ET.offset_x = dx * params.pix_x;
                        ET.offset_y  = dy * params.pix_y;
                        % offset = sqrt(offset_x^2 + offset_y^2);
                        ET = CheckFixation_VTD_egeVersion(ET, trial, params, output, center);
                    end
                end

            end
            
            
            if ET.TrialAborted == true
                flag = 1;
                if params.EEG == 1 
                     SendSignal(params.ioObj, params.address,55); %trial aborted trigger
                     triggerout = 1;
                     triggertime = GetSecs;
                else
                end
                %% Display April tags on PTB screen to detect surface via PupilLabs
                Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                Screen('DrawText', params.window, 'FIXATE PLEASE!', (params.screenX_pix/2)-100, (params.screenY_pix/2), params.black);
                Screen('Flip',params.window);
                break
            end
            
            
            % Timeout after state2
            if GetSecs>=output(trial).warningOnset+delayBeforeStimulus
                state = 3;
                break;
            else
            end
            
        end
        
        if flag == 1 || earlyResponseFlag==1
            fprintf('This trial is aborted at state # %d',state);
            if random_trl_order(trial) == 1
                stimPresentCount=stimPresentCount-1;
            elseif random_trl_order(trial) == 2
                stimAbsentCount=stimAbsentCount-1;
            end
            trial = trial - 1;
            continue
        end
        
  
        %% STATE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Send EEG trigger code for trial start
    
        
        if random_trl_order(trial) == 1 % For stimulus present trials
            % Apply the circular aperture to the
            % grating with the contrast recommended
            % by the QUEST algorithm
            params.gratingannulus = params.grey + params.grey * currentContrast * params.gratingtex;
            params.gratingannulus(params.annulusaperture == 0) = params.grey; % Set areas outside the annulus to gray
            
            % Creating the final texture from the grating with gray outside the annulus with PTB
            params.gratingannulustex = Screen('MakeTexture', params.window, params.gratingannulus);
            Screen('DrawTexture', params.window, params.gratingannulustex);
            Screen('FillOval', params.window, [0.5 0 0], params.centeredcircleRectOUT, params.maxDiameterOUT);
            Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
            Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, [0.5 0 0], [], 1);
            %Draw the peripheral white square for the photodiode
            Screen('DrawDots', params.window, [params.dotXpos1, params.dotYpos1], 10, [0.7 0.7 0.7], [], 4);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            % Get timestamp
            vbl = Screen('Flip', params.window);
            output(trial).stimulusOnset = vbl;
            if params.EEG == 1  %we only need this in stimulus present trials
                SendSignal(params.ioObj, params.address, 3); %3 = stimulus onset
                triggerout = 1;
                triggertime = GetSecs;
            else
            end
                            
            if params.eyetracking == 1 % Send Annotation from Matlab to Pupil Labs
                if params.whicheyetracker == 2
                    zmq.core.send(socket, uint8('t'));
                    currentTime_bytes = zmq.core.recv(socket,20);
                    % fprintf('first measure: %s\n', char(currentTime_bytes));
                    currentTimeChar = char(currentTime_bytes);
                    currentTimeNum = str2num(currentTimeChar);
                    keys_start = {'topic', 'label', 'timestamp', 'duration'};
                    values_start = {'annotation.StimulusOnset', strcat('Trial #', num2str(trial),' Stimulus Onset'), currentTimeNum, 1.0};
                    start_annotation = containers.Map(keys_start, values_start);
                    send_annotation(socket, start_annotation);
                    result = zmq.core.recv(socket);
                end
            end      
        elseif random_trl_order(trial) == 2 % For stimulus present trials
            
            % Draw stimulus and wait for stimulus presentation
            % duration but no stimulus is actually
            % presented
            % Screen('DrawTexture', params.window, params.gratingannulustex);
            % Get stimulus onset time
            output(trial).stimulusOnset = GetSecs;
            
        end
        
        flag4=0;
        earlyResponseFlag = 0;
        while state == 3
            [pressed, firstPress] = KbQueueCheck;
            
            if pressed
                % Extract information about the key presses
                [~,pressedKeys] = find(firstPress);
                
                % Check if any of the response keys are pressed
                if any(ismember(pressedKeys, KbName(responseKeys)))
                    % Identify the pressed key
                    % idx = find(ismember(KbName(responseKeys), pressedKeys));
                    % response = responseKeys{idx}; % Assign the pressed key as response
                    Screen('TextSize',params.window, 30)
                    Screen('DrawText',params.window, 'WAIT QUESTION MARK TO RESPOND!',(params.screenX_pix/2)-100, params.screenY_pix*0.5, params.black);
%                     Screen('DrawText',params.window, 'Please wait till the question mark sign appears on screen to make your report.',(params.screenX_pix/2)-100, params.screenY_pix*0.55, params.black);
                    %% Display April tags on PTB screen to detect surface via PupilLabs
                    Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                    Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                    Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                    Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                    Screen('Flip',params.window);
                    WaitSecs(2);                    
                    flag4 = 1;

                end
            end
            if flag4 == 1
                earlyResponseFlag=1;
                fprintf('early response');
                break
            end
            
            % Timeout after [params.LCandHC] seconds
            if GetSecs>=output(trial).stimulusOnset+ stimulusDuration
                state = 4;
                break;
            else
            end
%             
        end
%         


        %% STATE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        state = 4;
        if random_trl_order(trial) == 1
            % Remove stimulus and wait for variable
            % delay duration after stimulus, during
            % this period the red fixation mark still
            % remains
            Screen('FillOval', params.window, [0.5 0 0], params.centeredcircleRectOUT, params.maxDiameterOUT);
            Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
            Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, [0.5 0 0], [], 1);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            Screen('Flip', params.window, output(trial).stimulusOnset + (stimulusFrames - 0.5) * params.flipinterval);
            output(trial).delayOnset = GetSecs;
            % Get delay onset time

        elseif random_trl_order(trial) == 2 % For stimulus present trials
            % Remove stimulus and wait for variable
            % delay duration after stimulus, during
            % this period the red fixation mark still
            % remains
            Screen('FillOval', params.window, [0.5 0 0], params.centeredcircleRectOUT, params.maxDiameterOUT);
            Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
            Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, [0.5 0 0], [], 1);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            Screen('Flip', params.window, output(trial).stimulusOnset + (stimulusFrames - 0.5) * params.flipinterval);
            % Get delay onset time
            output(trial).delayOnset = GetSecs;
        end

        flag5=0;
        earlyResponseFlag = 0;
        while state == 4
            [pressed, firstPress] = KbQueueCheck;
            
            if pressed
                % Extract information about the key presses
                [~,pressedKeys] = find(firstPress);
                
                % Check if any of the response keys are pressed
                if any(ismember(pressedKeys, KbName(responseKeys)))
                    % Identify the pressed key
                    % idx = find(ismember(KbName(responseKeys), pressedKeys));
                    % response = responseKeys{idx}; % Assign the pressed key as response
                    Screen('TextSize',params.window, 30)
                    Screen('DrawText',params.window, 'WAIT QUESTION MARK TO RESPOND!',(params.screenX_pix/2)-100, params.screenY_pix*0.5, params.black);
%                     Screen('DrawText',params.window, 'Please wait till the question mark sign appears on screen to make your report.',(params.screenX_pix/2)-100, params.screenY_pix*0.55, params.black);
                    %% Display April tags on PTB screen to detect surface via PupilLabs
                    Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                    Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                    Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                    Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                    Screen('Flip',params.window);
                    WaitSecs(2);                    
                    flag5 = 1;

                end
            end
            if flag5 == 1
                earlyResponseFlag=1;
                fprintf('early response');
                break
            end
            % If eyetracking
            if params.eyetracking == 1
                if params.whicheyetracker == 2
                    if params.checkfix == 1
                        cursample = cursample + 1;
                        % Get gaze positions from pupil labs
                        % (filtermessages)
                        [topic, note] = recv_message(sub_socket, 2500);
                        if ~isequal(note, false)  % test for valid message
                            gaze_positions{1,iter}= [topic, note('gaze_on_surfaces')];  % print pupil norm_pos
                            for j=2:size(gaze_positions{1,iter},2)
                                gazecoordinates = gaze_positions{1,iter}{1,j}('norm_pos');
                                x(j) = gazecoordinates{1,1};
                                y(j) = gazecoordinates{1,2};
                                
                            end
                            %This alternative is decreasing the temporal...
                            %...resolution in case matlab fails.
                        end
                        
                        ET.CurrentGaze(1,1) = x(j);
                        ET.CurrentGaze(1,2) = y(j);
                        gazes(end+1,1) = ET.CurrentGaze(1,1);
                        gazes(end,2) = ET.CurrentGaze(1,2);
                        ET.currentTime = GetSecs();
                        ET.eyetime(cursample,1) = ET.currentTime;
                        dx = abs(x(j) - center(trial,1));
                        dy = abs(y(j) - center(trial,2));

                        ET.offset_x = dx * params.pix_x;
                        ET.offset_y  = dy * params.pix_y;
                        % offset = sqrt(offset_x^2 + offset_y^2);
                        ET = CheckFixation_VTD_egeVersion(ET, trial, params, output, center);
                    end
                end
            end
            
            if ET.TrialAborted == true
                flag = 1;
                if params.EEG == 1 
                     SendSignal(params.ioObj, params.address,55); %trial aborted trigger
                     triggerout = 1;
                     triggertime = GetSecs;
                else
                end
                %% Display April tags on PTB screen to detect surface via PupilLabs
                Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
                Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
                Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
                Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
                Screen('DrawText', params.window, 'FIXATE PLEASE!', (params.screenX_pix/2)-100, (params.screenY_pix/2), params.black);
                Screen('Flip',params.window);
                break
            end
            
            % EEG trigger off
            %             if triggerout == 1 && GetSecs >= triggertime + params.triggerofftime
            %                 SendSignal(0);
            %                 triggerout = 0;
            %                 triggertime = NaN;
            %             else
            %             end
            
            % Timeout after [params.LCandHC] seconds
            if GetSecs>=output(trial).delayOnset+delayAfterStimulus
                state = 5;
                break;
            else
            end
            
        end
        
        if earlyResponseFlag==1
            if params.EEG == 1 
                 SendSignal(params.ioObj, params.address,55); %trial aborted trigger
                 triggerout = 1;
                 triggertime = GetSecs;
            else
            end
        end
        
        if earlyResponseFlag==1 || flag==1
            fprintf('this trial is aborted at state # %d',state);
            if random_trl_order(trial) == 1
                stimPresentCount=stimPresentCount-1;
            elseif random_trl_order(trial) == 2
                stimAbsentCount=stimAbsentCount-1;
            end
            trial = trial - 1;
            continue
        end
        %% STATE 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Set the background to gray and display message
        % for participant to give their report and
        % confidence
        % Screen('FillRect', params.window, params.grey);
        % Screen('Flip', params.window);
        Report = '?';
        Screen('FillOval', params.window, params.black, params.centeredcircleRectOUT, params.maxDiameterOUT);
        Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);
        Screen('TextSize', params.window, 30);
        Screen('DrawText', params.window, Report, (params.screenX_pix/2)-10, (params.screenY_pix/2)-10, params.black);
        %% Display April tags on PTB screen to detect surface via PupilLabs
        Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
        Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
        Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
        Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));        
        Screen('Flip', params.window);
        % WaitSecs(time_for_response);
        
        response = [];
        flag6 = 0;
        flagResponse=1;
        responseStartTime = GetSecs(); % Record response start time
        while GetSecs() - responseStartTime <= time_for_response
            [pressed, firstPress] = KbQueueCheck;
            
            if pressed
                % Extract information about the key presses
                [~,pressedKeys] = find(firstPress);
                
                % Check if any of the response keys are pressed
                if any(ismember(pressedKeys, KbName(responseKeys)))
                    % Identify the pressed key
                    idx = find(ismember(KbName(responseKeys), pressedKeys));
                    response = responseKeys{idx}; % Assign the pressed key as response
                    break; % Exit the loop upon response detection
                end

            end
        end
        
        if ~isempty(response)
            if  random_trl_order(trial) == 1
                if strcmp(response, keyHit)
                    % Update QUEST for 'hit'
                    outputSP(stimPresentCount).response = 1;
                    output(trial).response = 1;
                    hitCount = hitCount + 1;
                    
                    if params.EEG == 1  %we only need this in stimulus present trials
                        SendSignal(params.ioObj, params.address, 4); %4 = response 'hit'
                        triggerout = 1;
                        triggertime = GetSecs;
                    else
                    end
                    
                    if stimPresentCount <= 10 | strcmp (params.runtype, 'TR') %if we are in the test run, we want contrast update after every trial
                        q = QuestUpdate(q, currentContrast, 1);
                    else
                        last10 = [outputSP(stimPresentCount-10:stimPresentCount).response];
                        sumlast10=sum(last10);
                        if sumlast10/10<0.5
                            q = QuestUpdate(q, currentContrast, 1);
                        elseif sumlast10/10>0.8
                            q = QuestUpdate(q, currentContrast, 1);
                        else                           
                        end
                    end
                    
                elseif strcmp(response, keyMiss)
                    % Update QUEST for 'miss'
                    outputSP(stimPresentCount).response = 0;
                    output(trial).response = 0;
                    missCount = missCount + 1;
                    
                    if params.EEG == 1  %we only need this in stimulus present trials
                        SendSignal(params.ioObj, params.address, 5); %5 = response 'miss'
                        triggerout = 1;
                        triggertime = GetSecs;
                    else
                    end
                    
                    if stimPresentCount <= 10
                        q = QuestUpdate(q, currentContrast, 0);
                    else
                        last10 = [outputSP(stimPresentCount-10:stimPresentCount).response]
                        sumlast10=sum(last10)
                        if sumlast10/10<0.5
                            q = QuestUpdate(q, currentContrast, 0);
                        elseif sumlast10/10>0.8
                            q = QuestUpdate(q, currentContrast, 0);
                        else
                        end
                    end
                    
                end
            elseif random_trl_order(trial) == 2
                if strcmp(response, keyHit)
                    faCount = faCount + 1;
                    output(trial).response = 1;
                elseif strcmp(response, keyMiss)
                    output(trial).response = 0;   
                end
            end
        elseif isempty(response)
            flag6 = 1;
            Screen('TextSize',params.window, 30)
            Screen('DrawText',params.window, 'RESPOND ON TIME!',(params.screenX_pix/2)-100, params.screenY_pix*0.5, params.black);
            %% Display April tags on PTB screen to detect surface via PupilLabs
            Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
            Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
            Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
            Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
            Screen('Flip',params.window);
            WaitSecs(2);
        end

     
        if flag6 == 1
            flagResponse = 0;
        end

        if flagResponse == 0
            if params.EEG == 1 
                 SendSignal(params.ioObj, params.address,55); %trial aborted trigger
                 triggerout = 1;
                 triggertime = GetSecs;
            else
            end
            fprintf('this trial is aborted at state # %d',state);
            if random_trl_order(trial) == 1
                stimPresentCount=stimPresentCount-1;
            elseif random_trl_order(trial) == 2
                stimAbsentCount=stimAbsentCount-1;
            end
            trial = trial - 1;
            continue
        end

        WaitSecs(0.5);
        [pos_conf] = slideScale(params.window, 'Rate your confidence', params.resolution, {'LOW', 'AVERAGE', 'HIGH'},params, 'device','keyboard');
        %% Display April tags on PTB screen to detect surface via PupilLabs
        Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
        Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
        Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
        Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
        Screen('Flip',params.window);
        WaitSecs(0.5);
        output(trial).confidenceAnswer = pos_conf; %deviation from the center

        hitRates(trial) = (hitCount / stimPresentCount)*100; % Calculate hit rate
        falseAlarmRates(trial) = (faCount / stimAbsentCount)*100; % Calculate false alarm rate
        output(trial).accuracy = hitRates(trial);
        output(trial).falsealarm = falseAlarmRates(trial);
        fprintf('Current hit rate: %d',hitRates(trial));
        fprintf('\nCurrent false alarm rate: %d',falseAlarmRates(trial));
        
%         % % Update data arrays
%         x1 = [x1, trial];   % X-axis as index
%         y1 = [y1, hitRates(trial)];
%         x2 = [x2, trial];   % X-axis as index
%         y2 = [y2, falseAlarmRates(trial)];
%         x3 = [x3, trial];
%         y3 = [y3, currentContrast];
%         
%         % Plot online for the first figure
%         figure(1);  % Activate the first figure
%         subplot(2, 1, 1);
%         plot(x1, y1, 'r-', 'LineWidth', 2);  % Plot updated data for plot 1
%         hold on;  % Hold the current plot to add more data
%         plot(x1, y2, 'b--', 'LineWidth', 2);  % Plot y2 in red dashed line
%         hold off
%          % Release the current plot
%         xlabel('Trial');
%         ylabel('Hits and False Alarms');
%         ylim([0,1]);
%         title('QUEST modulated accuracy');
%         drawnow;  % Refresh the plot
% 
%         
%         % Plot online for the second figure
%         % figure(2);  % Activate the second figure
%         subplot(2, 1, 2);
%         plot(x3, y3, 'k-', 'LineWidth', 2);  % Plot updated data for plot 2
%         xlabel('Trial');
%         ylabel('Grating Contrast');
%         ylim([0,0.3]);
%         title('QUEST modulated contrast');
%         drawnow;  % Refresh the plot
        
        
        
        
        % EEG trigger off
%         if triggerout == 1 && GetSecs >= triggertime + params.triggerofftime
%             SendSignal(0);
%             triggerout = 0;
%             triggertime = NaN;
%         else
%         end;
        
        
        fprintf('\nThis trial is successfully completed, trial number: %d',trial);
        
        if params.checkfix == 1
            newMean = [mean(gazes(:,1)) mean(gazes(:,2))]; %after a successful completion of a trial, you get the mean gaze locations to re-define it as the center!
            output(trial).meanGazes = newMean;
        end
        
        if trial == params.numtrials
            break
        end

    end; % end of trial loop

    WaitSecs(2);  
    KbQueueStop;	% Stop delivering events to the queue

    [pressed, firstPress] = KbQueueCheck;
       
% [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck;

    % Create key presses struct
    for i = 1:KbEventAvail
        [evt, n] = KbEventGet;
        if i == 1
            keyboard_events = evt;
        else
            keyboard_events(end+1) = evt;
        end;
    end;


end % end of blocks loop
%% Close properly

% Show end message for the subject
Screen('DrawText', params.window, 'Done - Thank you!', (params.screenX_pix/2)-100, (params.screenY_pix/2), params.black);
Screen('Flip', params.window);
WaitSecs(2); 

% Show mouse cursor again
ShowCursor;

% % Calculate hit rate and false alarm rate as percentages
overallhitRate = (hitCount / stimPresentCount) * 100;
overallfalseAlarmRate = (faCount / stimAbsentCount) * 100;


% Save data to .mat file
cd(params.savepath);

save(params.filename,'params','output','keyboard_events','ET','gazes','-v7.3');

disp('Data saved');

KbQueueRelease;
ListenChar(0);

if strcmp (params.runtype, 'TR')| strcmp (params.runtype, 'TC') 
% Display the estimated threshold
    disp(['Estimated Threshold from the Calibration Run: ' num2str(QuestMean(q))]);
elseif strcmp (params.runtype, 'R')| strcmp (params.runtype, 'H') 
    disp(['Estimated Threshold from the Last Run: ' num2str(QuestMean(q))]);
end

% Close Eyelink eyetracker
if params.eyetracking == 1
   if params.whicheyetracker == 1 
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        % download data file
             try
                fprintf('Receiving data file ''%s''\n', params.eyedatafile);
                status=Eyelink('ReceiveFile');
                if status > 0
                    fprintf('ReceiveFile status %d\n', status);
                end
                if 2==exist(params.eyedatafile, 'file')
                    fprintf('Data file ''%s'' can be found in ''%s''\n', params.eyedatafile, pwd );
                end
            catch rdf
                fprintf('Problem receiving data file ''%s''\n', params.eyedatafile);
                rdf;
            end

        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        Eyelink('shutdown')
   elseif params.whicheyetracker == 2
        %stop the recording of Pupil Labs
        zmq.core.send(socket, uint8('r'));
        result = zmq.core.recv(socket);
        fprintf('Recording stopped: %s\n', char(result));
        
        % disconnect sub socket
        zmq.core.disconnect(sub_socket, sub_endpoint);
        zmq.core.close(sub_socket);
        fprintf('Disconnected from SUB: %s\n', sub_endpoint);
        
        zmq.core.disconnect(socket, endpoint);
        zmq.core.close(socket);

%         zmq.core.ctx_shutdown(ctx);
%         zmq.core.ctx_term(ctx);
   end 
end



% Find indices of stimulus-absent trials
stimulusAbsentIdx = find(random_trl_order == 2);
% Generate plots
% Plot for Accuracy (hit rate and false alarm rate vs Trial Number)
trialNumbers = 1:length(random_trl_order);
figure;
subplot(2, 1, 1);
plot(trialNumbers, hitRates, 'Color', [0.6, 0, 0], 'LineWidth', 2); % Plot hit rate
hold on;
plot(trialNumbers, falseAlarmRates, 'Color', [0, 0.6, 0], 'LineWidth', 2); % Plot false alarm rate
xlabel('Trial Number');
ylabel('Accuracy (%)');
% legend('Hit Rate', 'False Alarm Rate', 'Location','eastoutside');
title('Accuracy');
xlim([1, length(trialNumbers)]);

% Plot the stimulus-absent trials individually with a grey background
for i = 1:length(stimulusAbsentIdx)
    idx = stimulusAbsentIdx(i);
    
    % Check if this is not the last stimulus-absent trial
    if i < length(stimulusAbsentIdx)
        % Define x-coordinates for the grey background area
        x = [idx, trialNumbers(idx+1), trialNumbers(idx+1), idx];
    else
        % For the last stimulus-absent trial
        if idx == trialNumbers(end)
            x = [idx, trialNumbers(end), trialNumbers(end), idx];
        else
            x = [idx, trialNumbers(idx+1), trialNumbers(idx+1), idx];
        end
    end
    
    % Define y-coordinates for the grey background area (0 to 0.25 for contrast)
    y = [0, 0, 100, 100];
    
    % Plot the grey background area for this stimulus-absent trial
    patch(x, y, 'k', 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);
end


% Plot for currentContrast vs Trial Number
subplot(2, 1, 2);
plot(trialNumbers, contrasts, 'Color', [0.6, 0, 0], 'LineWidth', 2);
xlabel('Trial Number');
ylabel('Stimulus Contrast');
title('QUEST recommended contrast');
ylim([0, 0.25]); % Set y-axis limits for contrast plot
xlim([1, length(trialNumbers)]);

% Plot the stimulus-absent trials individually with a grey background
for i = 1:length(stimulusAbsentIdx)
    idx = stimulusAbsentIdx(i);
    
    % Check if this is not the last stimulus-absent trial
    if i < length(stimulusAbsentIdx)
        % Define x-coordinates for the grey background area
        x = [idx, trialNumbers(idx+1), trialNumbers(idx+1), idx];
    else
        % For the last stimulus-absent trial
        if idx == trialNumbers(end)
            x = [idx, trialNumbers(end), trialNumbers(end), idx];
        else
            x = [idx, trialNumbers(idx+1), trialNumbers(idx+1), idx];
        end
    end
    
    % Define y-coordinates for the grey background area (0 to 0.25 for contrast)
    y = [0, 0, 0.25, 0.25];
    
    % Plot the grey background area for this stimulus-absent trial
    patch(x, y, 'k', 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);
end
hold off;

% Save the plots 
saveas(gcf, fullfile(params.savepathplot, params.filenameplot));

end

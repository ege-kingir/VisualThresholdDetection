%% VTD Main: With Cycling, or without a secondary task!
% Sukanya C, Ege Kingir (MSc Thesis, IMPRS)
% 12 Nov 2023

% SESSION OUTLINE: 
% 1. Forms, general instructions

% 2. Record 10 min baseline and measure blood pressure
% 3. Task Instructions; Do 6-8 blocks (87% stimulus present trials in each block): 
% With cycling, a total of 7 runs : 2 with low-resistance cycling (L), 2 with high-resistance cycling (H), 3 rest (R).
% Allow for enough break time between runs (esp. to allow heart rate to
% come back down to normal).
% In pseudo-random order: LRHRLRH or HRLRHRL.


% 4. Record 8 min baseline and measure blood pressure


%% Hardware:
% 1) stimulus presentation computer (MATLAB-psychtoolbox)
% 2) physiological recording computer (BrainVision Recorder)
% 3) Pupil Labs eye-tracking computer, and the
% eye-tracking glasses (Pupil Core)
% 4) Trigger streamline: Ethernet between computer #1 and #3, TTL (USB-2
% adapter BrainVision) between computer #1 and #2


%% CHANGE LOG
% Following changes made by XX on XX.XX.XXXX
%%%


%% FUNCTIONS CALLED IN THIS SCRIPT / RELEVANT FILES
% params_VTD.m
% startupEyetracking.m
% SendSignal.m
% display_cycle_clock.m (only for blocks with cycling)
% trials_VTD.m


%%
clear all;
close all;
clc; 

par.screenNumber = input('Screen number (test using ScreenTest before starting experiment) '); % almost always 1

try
%% USER INPUT
fprintf('VTD \n\nRun type options: \nHigh-resistance Cycling 1 (H) \nLow-resistance Cycling (L) \nRest (R) \nTest with cycling (TC) \nTest with rest (TR) \n\n\n');


% Run Sequence A to F or test
par.runtype = input('Run type (e.g. TR = Test with Rest, H = High-resistance cycling, TC = Test with Cycling) ', 's');

if strcmp(par.runtype, 'R') | strcmp(par.runtype, 'L') | strcmp(par.runtype, 'H')
    % Anonymous subject number (e.g. 001, 020, 047, 194)
    par.subject = input('Subject (e.g. 001, 020, 047, 194) ', 's');
    
    % Number of the current run (e.g. 1-7)
    par.runnumber = input('Run number (e.g. between 1-7) ');
%     
    % Eyetracking using EyeLink eyetracker on or off (binary)
    par.eyetracking = input('Last question: Eyetracking (1 = on, 0 = off) ');
    
%     if par.eyetracking == 1
        par.whicheyetracker = input('Which eyetracker do you want to use? (Eyelink = 1, Pupillabs = 2) ');
        par.checkfix = input('Check fixation breaks? (1 = on, 0 = off) ');
%     end;
    par.EEG = input('0=off, 1=on' ); % 1 = on
    
    
elseif strcmp(par.runtype, 'TR') | strcmp(par.runtype, 'TC')
    % Anonymous subject number (e.g. 001, 020, 047, 194)
    par.subject = input('Subject (e.g. 001, 020, 047, 194) ', 's');
    
    % Number of the current run (e.g. 1-7)
    par.runnumber = input('Run number (e.g. between 1-3) ');
%     
    % Eyetracking using EyeLink eyetracker on or off (binary)
    par.eyetracking = input('Last question: Eyetracking (1 = on, 0 = off) ');
    

    par.whicheyetracker = input('Which eyetracker do you want to use? (Eyelink = 1, Pupillabs = 2) ');
    par.checkfix = input('Check fixation breaks? (1 = on, 0 = off) ');

    
end

%% RUN EXPERIMENT

% Get params
par = params_VTD_egeVersion(par);

% Main trial function
[out, key_events, par, EyeTracking, gaze] = trials_VTDnew_egeVersion(par);

% Close open PTB windows
Screen('CloseAll')
return;

catch ERR
   % Show mouse cursor again
    ShowCursor;
    
    KbQueueRelease;
    ListenChar(0);
    
    rethrow(ERR) % where ERR is the variable name under which the error is stored
    
    % Close open PTB windows
    Screen('CloseAll')
    return;
end;

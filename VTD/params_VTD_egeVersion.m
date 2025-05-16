function [params] = params_VTD_egeVersion(params)
% Parameters function for MAIN_TVD.m

% Sukanya C, Ege Kingir (MSc Thesis, IMPRS)
% 07 November 2023


%% CHANGE LOG
% Following changes made by XX on XX.X.XXXX
% 30.11.2023 EK changes:
% 1) 
%%


%% FUNCTIONS CALLED IN THIS SCRIPT / RELEVANT FILES


%% Psychtoolbox parameters

% Here we call some default settings for setting up PTB
addpath('C:\toolbox\Psychtoolbox\PsychOneliners\')
PsychDefaultSetup(2);

% Screen('Preference', 'SkipSyncTests', 1); % Skip sync tests for testing purposes

params.black = BlackIndex(params.screenNumber);
params.white = WhiteIndex(params.screenNumber);
% Do a simply calculation to calculate the luminance value for grey. This
% will be half the luminace value for white
params.grey = params.white / 2;

% Open grey window
[params.window, params.resolution] = PsychImaging('OpenWindow', params.screenNumber, params.grey);

% Enable alpha blending for anti-aliasing (smooth edges / lines)
Screen('BlendFunction', params.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Size of the on screen window in pixels
[params.screenX_pix, params.screenY_pix] = Screen('WindowSize', params.window);

% Centre of screen in pixels
[params.screenXcentre_pix, params.screenYcentre_pix] = RectCenter(params.resolution);
params.centrecoords_pix = [params.screenXcentre_pix params.screenYcentre_pix];

Screen('TextFont',params.window,'Arial');
% Screen('TextSize',params.window, 30);

% Frame duration / flip duration
params.flipinterval = Screen('GetFlipInterval', params.window);

% Screen refresh rate in Hz
params.refreshrate = round(1/params.flipinterval); % from frame duration
params.hz = Screen('NominalFrameRate', params.window); % actual frame rate

% Maximum priority level
topPriorityLevel = MaxPriority(params.window);
Priority(topPriorityLevel);

% Flip to black screen
Screen('Flip', params.window);


%%
%
% Seed the random number generator. Here we use an older way to be compatible with older systems.
rng('shuffle');


%% Parameters of physical setup

% Distance to screen in cm
params.dist2screen_cm = 70;

% Screen params.screenxres in cm
params.screenxsize_cm = 30;

% Screen params.screenyres in cm
params.screenysize_cm = 20;

% X and y resolution of the screen in pixel
params.screenxres_pix = params.resolution(3);
params.screenyres_pix = params.resolution(4);

% Calculate number of pixels per degree
params.pixperdeg = round((params.dist2screen_cm/params.screenxsize_cm) * ...
    (sin((pi * 1)/180) / sin((pi * (90 - 1))/180)) * params.screenxres_pix);

% Center of the screen in pixel
params.xcenter_pix = params.resolution(3)/2;
params.ycenter_pix = params.resolution(4)/2;

%% EEG triggers

% Time after trigger signal trigger can be reset to zero (seconds)
% params.triggerofftime = 0.1;

% Triggers key
% 1 = SP trial onset
% 2 = SA trial onset
% 3 = Stimulus onset


% Trigger object and address params
params.ioObj = io64;
params.address = hex2dec('DFB8');

%% File saving parameters
%
% Folder in which the data are saved
params.savepath = 'C:\Users\sheumue\Desktop\Sukanya\VTD\Data';
params.savepathplot = 'C:\Users\sheumue\Desktop\Sukanya\VTD\Plots';
% Current experiment date and time
params.exptime = datestr(now);

params.filename = sprintf('VTD_Subj%s_Run%d_Type%s_%s.mat',params.subject, params.runnumber, params.runtype, datestr(params.exptime,30));
params.filenameplot = sprintf('VTD_Subj%s_Run%d_Type%s_%s.png',params.subject, params.runnumber, params.runtype, datestr(params.exptime,30));


%% Trial numbers, order, and timings

% An experimental session consisted of six to eight blocks with 80
% stimulus-present trials and 12 stimulus-absent trials interleaved 
% randomly per block. 


% % Total number of trials in a run
if strcmp (params.runtype, 'TR')
    params.numblocks = 1; % Random number of blocks between 6 and 8, each block needs to be run separately
    params.SPnumtrials = 30;
    params.SAnumtrials = 5;
    params.numtrials = params.SPnumtrials + params.SAnumtrials; % total number of trials in a run
elseif strcmp (params.runtype, 'R') | strcmp (params.runtype, 'L') | strcmp (params.runtype, 'H')
    % % Conditions: 1 = Stimulus Present (SP), 2 = Stimulus Absent (SA)
    params.numblocks = 1; % Random number of blocks between 6 and 8, each block needs to be run separately
    params.SPnumtrials = 61;
    params.SAnumtrials = 9;
    params.numtrials = params.SPnumtrials + params.SAnumtrials; % total number of trials in a run
    params.cycledur = 60;
elseif strcmp (params.runtype, 'TC')
    params.numblocks = 1; % Random number of blocks between 6 and 8, each block needs to be run separately
    params.SPnumtrials = 30;
    params.SAnumtrials = 5;
    params.numtrials = params.SPnumtrials + params.SAnumtrials; % total number of trials in a run
    params.cycledur = 60;
end
% Randomly interleave stimulus-present and stimulus-absent trials
params.stimuluspresentProb = params.SPnumtrials/params.numtrials; % stimulus present probability is 87 % in the total trials
% params.trialorder = Shuffle([ones(1, params.numtrials * params.stimuluspresentProb), ...
%                               repmat(2, 1, params.numtrials * (1 - params.stimuluspresentProb))]);
params.trialorder = Shuffle([ones(1, params.SPnumtrials), ...
                              repmat(2, 1, params.SAnumtrials)]);

% Range for fixation duration in seconds
params.fixationDurRange = [0.5, 0.7];
params.fixdur = 1;
% Range for variable delay before stimulus presentation
params.vardelayDurRange = [1.5, 1.7];
% Stimulus duration in seconds
params.stimuluspresentDur = round(0.05*params.refreshrate);
% Range for intertrial interval in seconds
params.intertrialDurRange = [1.5, 2];
     
%% Stimulus size and location parameters

% Stimuli consisted of a grating annulus (spatial frequency, five cycles per degree of visual angle;
% inner and outer radius, 2.5° and 3° of visual angle, respectively; orientation chosen among 20
% equally spaced between 0° and 180°, with cardinal orientations (±30° around the horizontal and
% vertical orientations) being excluded). The fixation mark was a black dot (radius, 0.15° of visual angle)
% surrounded by a black circle (radius, 0.5° of visual angle) at the center of the screen. All stimuli were
% presented on a gray background (luminance, 4.6 cd m−2) at a viewing distance of 0.8 m.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GRATING ANNULUS                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Following parameters in degrees of visual angles %%%

params.gratingFrequency = 5; % Spatial frequency (cycles per degree; Park et al., 2014)
params.innerRadiusDeg = 3; % Inner radius of the annulus (degrees) default is 2.5
params.outerRadiusDeg = 3.6; % Outer radius of the annulus (degrees) default is 3
params.contrast = 1; % Contrast of the grating, also the parameter subject to modulation to find detection threshold

%%% Corresponding pixel values for stimulus parameters %%%
params.gratingFrequency_cycperpix = 5/params.pixperdeg;
params.innerRadius_pix = params.innerRadiusDeg*params.pixperdeg;
params.outerRadius_pix = params.outerRadiusDeg*params.pixperdeg;

%%% Determining the random orientation for grating lines %%%
params.gratingSize = params.resolution(3); % Use the entire screen dimensions

% Define the number of lines you want (i.e., 20 lines; default from Park et al., 2014)
params.gratingnumLines = 20;

% Define the excluded angle ranges (±30° around the horizontal and
% vertical orientations excluded; default from Park et al., 2014)
params.gratingexcludedRanges = [0, 30; 60, 120; 150, 180];

% Initialize the angle variable
params.gratingangle = [];

        % Generate an angle that is not in the excluded ranges
        while isempty(params.gratingangle)
            % Generate a set of equally spaced angles between 0 and 180 degrees
            params.angles = linspace(0, 180, params.gratingnumLines);
        
            % Randomly select one angle from the set
            randInd = randperm(length(params.angles), 1);
            params.gratingangle = params.angles(randInd);
        
            % Check if the angle falls within any excluded range, then set
            % it to empty and try to generate a different angle
            for angle_idx = 1:size(params.gratingexcludedRanges, 1)
                if params.gratingangle >= params.gratingexcludedRanges(angle_idx, 1) && params.gratingangle <= params.gratingexcludedRanges(angle_idx, 2)
                    params.gratingangle = [];
                    break; % Regenerate angle
                end
            end
        end

%%% Preparing the grating and restricting it to an annulus shape %%%
% Create a grating texture
[x, y] = meshgrid(-params.gratingSize:params.gratingSize, -params.gratingSize:params.gratingSize);
params.gratingtex = cosd(360 * params.gratingFrequency_cycperpix * (cosd(params.gratingangle) * x + sind(params.gratingangle) * y));

% Create a circular aperture for the annulus
[x, y] = meshgrid(-params.gratingSize:params.gratingSize, -params.gratingSize:params.gratingSize);
params.annulusaperture = double((x.^2 + y.^2 >= params.innerRadius_pix^2) & (x.^2 + y.^2 <= params.outerRadius_pix^2));

% % Apply the circular aperture to the grating
% params.gratingannulus = params.grey + params.grey * params.contrast * params.gratingtex;
% params.gratingannulus(params.annulusaperture == 0) = params.grey; % Set areas outside the annulus to gray

% Creating the final texture from the grating with gray
% outside the annulus with PTB
% params.gratingannulustex = Screen('MakeTexture', params.window, params.gratingannulus);

%% params for the white dot for the photo diode
params.dotXpos1 = 0.5*params.screenX_pix;
params.dotYpos1 = 0.03*params.screenY_pix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   FIXATION CIRCLE                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Following parameters in degrees of visual angles %%%
params.fixationSizeDeg = 0.15; % Size of the fixation dot (degrees)
params.fixationCircleSizeDegOUT = 0.42; % Size of the outer fixation circle (degrees)
params.fixationCircleSizeDegIN = 0.35; % Size of the inner fixation circle (degrees)

%%% Corresponding pixel values for stimulus parameters %%%
params.fixationSize_pix = params.fixationSizeDeg*params.pixperdeg;
params.fixationCircleSizeOUT_pix = params.fixationCircleSizeDegOUT*params.pixperdeg;
params.fixationCircleSizeIN_pix = params.fixationCircleSizeDegIN*params.pixperdeg;

%%% Creating the outer circle for fixation mark %%%
params.circleRectOUT = [params.xcenter_pix - params.fixationCircleSizeOUT_pix, params.ycenter_pix - params.fixationCircleSizeOUT_pix, params.xcenter_pix + params.fixationCircleSizeOUT_pix, params.ycenter_pix + params.fixationCircleSizeOUT_pix];
% For Ovals we set a miximum diameter up to which it is perfect for
params.maxDiameterOUT = max(params.circleRectOUT) * 1.01;
% Center the rectangle on the centre of the screen
params.centeredcircleRectOUT = CenterRectOnPointd(params.circleRectOUT, params.xcenter_pix, params.ycenter_pix);
%%% Screen('FillOval', params.window, params.black, params.centeredcircleRectOUT, params.maxDiameterOUT);

%%% Creating the inner circle for fixation mark %%%
params.circleRectIN = [params.xcenter_pix - params.fixationCircleSizeIN_pix, params.ycenter_pix - params.fixationCircleSizeIN_pix, params.xcenter_pix + params.fixationCircleSizeIN_pix, params.ycenter_pix + params.fixationCircleSizeIN_pix];
% For Ovals we set a miximum diameter up to which it is perfect for
params.maxDiameterIN = max(params.circleRectIN) * 1.01;
% Center the rectangle on the centre of the screen
params.centeredcircleRectIN = CenterRectOnPointd(params.circleRectIN, params.xcenter_pix, params.ycenter_pix);
%%% Screen('FillOval', params.window, params.grey, params.centeredcircleRectIN, params.maxDiameterIN);

% Drawing the central dot for fixation mark
%%% Screen('DrawDots', params.window, [params.xcenter_pix, params.ycenter_pix], params.fixationSize_pix, params.black, [], 1);


%% Color scheme

% Target stimulus color in RGB
params.annulusColor = [255 255 255]; % Default full-contrast white, needs to later be adjusted to threshold

% Fixation dot color in RGB
params.fixationColor = [0 0 0]; % Default black; Park et al., 2014

% Fixation circle color in RGB
params.fixationCircleColor = [0 0 0]; % Default black; Park et al., 2014

%% April tags
params.tag1 = imread('AprilTag tag36h11_1_highRes-1_lowBright_lowContrast.png');
params.tag2 = imread('AprilTag tag36h11_2_highRes-1_lowBright_lowContrast.png');
params.tag3 = imread('AprilTag tag36h11_3_highRes-1_lowBright_lowContrast.png');
params.tag4 = imread('AprilTag tag36h11_4_highRes-1_lowBright_lowContrast.png');

%make textures from images
params.txt1 = Screen('MakeTexture', params.window, params.tag1);
params.txt2 = Screen('MakeTexture', params.window, params.tag2);
params.txt3 = Screen('MakeTexture', params.window, params.tag3);
params.txt4 = Screen('MakeTexture', params.window, params.tag4);

% Get the sizes of the images
[params.s11, params.s21, params.s31] = size(params.tag1);
[params.s12, params.s22, params.s32] = size(params.tag2);
[params.s13, params.s23, params.s33] = size(params.tag3);
[params.s14, params.s24, params.s34] = size(params.tag4);

% Get the aspect ratio of the image. We need this to maintain the aspect
% ratio of the image when we draw it different sizes. Otherwise, if we
% don't match the aspect ratio the image will appear warped / stretched
params.aspectRatio1 = params.s21 / params.s11;
params.aspectRatio2 = params.s22 / params.s12;
params.aspectRatio3 = params.s23 / params.s13;
params.aspectRatio4 = params.s24 / params.s14;

params.imageHeights = params.screenY_pix/6;
params.imageWidths = round(params.imageHeights*params.aspectRatio1);

params.numImages =4;
params.pixelsOnTheSide = 100;
% params.x = [0.10 0.10 0.90 0.90];
% params.y = [0.15 0.85 0.15 0.85];
params.x = [round(params.imageWidths/2)+params.pixelsOnTheSide round(params.imageWidths/2)+params.pixelsOnTheSide params.screenX_pix-(round(params.imageWidths/2)+params.pixelsOnTheSide) params.screenX_pix-(round(params.imageWidths/2)+params.pixelsOnTheSide)];
params.y = [round(params.imageHeights/2)+params.pixelsOnTheSide params.screenY_pix-(round(params.imageHeights/2)+params.pixelsOnTheSide) round(params.imageHeights/2)+params.pixelsOnTheSide params.screenY_pix-(round(params.imageHeights/2)+params.pixelsOnTheSide)];
for i=1:params.numImages
    params.theRect = [0 0 params.imageWidths params.imageHeights];
    params.dstRects(:, i) = CenterRectOnPointd(params.theRect, params.x(i),...
    params.y(i));
end

params.pix_x = params.screenX_pix - 2*params.pixelsOnTheSide; % Number of horizontal pixels within pupillabs defined surface
params.pix_y = params.screenY_pix - 2*params.pixelsOnTheSide; % Number of vertical pixels within pupillabs defined surface

%% params for eye-tracking / fixation control
params.fixRadiusdeg = 2; % Window in degrees of visual angles for checking for fixation
params.fixRadius = params.fixRadiusdeg*params.pixperdeg; % Window in pixels for checking for fixation

%reliance on previous fixation centers vs. the most recent fixation
params.weightOld = 0.90;
params.weightNew = 0.10;

%% EEG triggers

% EEG triggers on or off (binary)
% params.EEG = 1; % 1 = on

% Time after trigger signal trigger can be reset to zero (s)
params.triggerofftime = 0.1;

% Triggers key
% 1 = SP
% 2 = SA
% 3 = Blank fixation
     
end
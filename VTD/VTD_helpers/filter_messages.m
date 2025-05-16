% pupil_remote_control.m

% (*)~----------------------------------------------------------------------------------
%  Pupil Helpers
%  Copyright (C) 2012-2016  Pupil Labs
% 
%  Distributed under the terms of the GNU Lesser General Public License (LGPL v3.0).
%  License details are in the file license.txt, distributed as part of this software.
% ----------------------------------------------------------------------------------~(*)

clear;clc
%% Open a window with April Tags
params.screenNumber = input('Screen number (test using ScreenTest before starting experiment) '); % almost always 1

PsychDefaultSetup(2);

Screen('Preference', 'SkipSyncTests', 1); % Skip sync tests for testing purposes
% % Define black and white
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


%% Pupil Labs calibration
        
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

        Screen('FillOval', params.window, [0.5 0.5 0.5],centeredRect1, max(baseRect1)*1.01);
        Screen('FillOval', params.window, [1 1 1],centeredRect2, max(baseRect2)*1.01);
        Screen('FillOval', params.window, [0 0 0],centeredRect3, max(baseRect3)*1.01);
        Screen('FillOval', params.window, [1 1 1],centeredRect4, max(baseRect4)*1.01);
        Screen('FillOval', params.window, [0 0 0],centeredRect5, max(baseRect5)*1.01);

        %fixation cross
        Screen('DrawLines', params.window, allCoords, lineWidthPix, [1 1 1], [params.screenXcentre_pix+x(stim) params.screenYcentre_pix+y(stim)],2);
        Screen('Flip', params.window);
        KbStrokeWait;

    end
    Screen('Flip', params.window);

end
WaitSecs(5);
Screen('Flip',params.window);

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
%% Display April tags on PTB screen to detect surface via PupilLabs
Screen('DrawTextures', params.window, params.txt1, [], params.dstRects(:,1));
Screen('DrawTextures', params.window, params.txt2, [], params.dstRects(:,2));
Screen('DrawTextures', params.window, params.txt3, [], params.dstRects(:,3));
Screen('DrawTextures', params.window, params.txt4, [], params.dstRects(:,4));
Screen('Flip', params.window);

% Setup zmq context and remote helper
ctx = zmq.core.ctx_new();
req_socket = zmq.core.socket(ctx, 'ZMQ_REQ');

% set timeout to 1000ms in order to not get stuck in a blocking
% mex-call if server is not reachable, see
% http://api.zeromq.org/4-0:zmq-setsockopt#toc19
zmq.core.setsockopt(req_socket, 'ZMQ_RCVTIMEO', 1000);

ip_address = '169.254.53.40';
req_port = '50020';
req_endpoint =  sprintf('tcp://%s:%s', ip_address, req_port);

fprintf('Connecting to REQ: %s\n', req_endpoint);
zmq.core.connect(req_socket, req_endpoint);

% Request sub port
zmq.core.send(req_socket, uint8('SUB_PORT'));
sub_port = char(zmq.core.recv(req_socket));
fprintf('Received sub port: %s\n', sub_port);

% Disconnect req socket
zmq.core.disconnect(req_socket, req_endpoint);
zmq.core.close(req_socket);
fprintf('Disconnected from REQ: %s\n', req_endpoint);

if isequal(sub_port, false)
    warning('No valid sub port received');
    return;  % exit script
end

% Create and connect sub socket
sub_endpoint =  sprintf('tcp://%s:%s', ip_address, sub_port);
sub_socket = zmq.core.socket(ctx, 'ZMQ_SUB');

% set timeout to 1000ms in order to not get stuck in a blocking
% mex-call if server is not reachable, see
% http://api.zeromq.org/4-0:zmq-setsockopt#toc19
zmq.core.setsockopt(sub_socket, 'ZMQ_RCVTIMEO', 1000);

fprintf('Connecting to SUB: %s\n', sub_endpoint);
zmq.core.connect(sub_socket, sub_endpoint);

% set subscriptions to topics
% recv just pupil/gaze/notifications
% zmq.core.setsockopt(sub_socket, 'ZMQ_SUBSCRIBE', 'pupil.');
% zmq.core.setsockopt(sub_socket, 'ZMQ_SUBSCRIBE', 'gaze');
zmq.core.setsockopt(sub_socket, 'ZMQ_SUBSCRIBE', 'surface');

% zmq.core.setsockopt(socket, 'ZMQ_SUBSCRIBE', 'notify.');
% zmq.core.setsockopt(socket, 'ZMQ_SUBSCRIBE', 'logging.');
% or everything else
% zmq.core.setsockopt(socket, 'ZMQ_SUBSCRIBE', '');
% [topic, note] = 
% surface_name = 'Surface 1';
figure
% for iter=1:100 % receive the first whatever messages
%     % messages that are longer than 1024 bytes will be tuncated
%     % and ignored. In this case note is set to false. Increase
%     % the buffer size if you experience this issue.
%     [topic, note] = recv_message(sub_socket, 2000);
%     
%     if ~isequal(note, false)  % test for valid message
%         gaze_positions {1,iter}= [topic, note('gaze_on_surfaces')];   % print pupil norm_pos
%         for j=2:size(gaze_positions{1,iter},2)
%             gazecoordinates = gaze_positions{1,iter}{1,j}('norm_pos');
%             x(j) = gazecoordinates{1,1};
%             y(j) = gazecoordinates{1,2};
%             plot(x(j),y(j),'p')
%             drawnow
%             hold on
%         end
%     end
% end

iter=1;
total_count=1;
for i=1:500
    [topic, note] = recv_message(sub_socket, 2500);
    if ~isequal(note, false)  % test for valid message
        gaze_positions{1,iter} = [topic, note('gaze_on_surfaces')];   % print pupil norm_pos
        for j=2:size(gaze_positions{1,iter},2)
            gazecoordinates = gaze_positions{1,iter}{1,j}('norm_pos');
            x(total_count) = gazecoordinates{1,1}; 
            y(total_count) = gazecoordinates{1,2};
%             plot(x(j),y(j),'p')
%             drawnow
%             hold on
            total_count = total_count+1;
        end
    end
    iter=iter+1;
%     clear gaze_positions gazecoordinates x y
end

% disconnect sub socket
zmq.core.disconnect(sub_socket, sub_endpoint);
zmq.core.close(sub_socket);
fprintf('Disconnected from SUB: %s\n', sub_endpoint);

zmq.core.ctx_shutdown(ctx);
zmq.core.ctx_term(ctx);

% [norm_pos_x(i,j),norm_pos_y(i,j)] = 
% for i = 1:1
%     for j = 2:size(gaze_positions{1,i},2)
%         confidence(i,j) = gaze_positions{1,i}{1,j}('confidence');
%         gaze = gaze_positions{1,i}{1,j}('norm_pos');
%         norm_pos_x(i,j) = gaze{1,1};
%         norm_pos_y(i,j) = gaze{1,2};
%     end
% end

for i=1:length(x)
    plot(x(i),y(i),'p');
    hold on
end
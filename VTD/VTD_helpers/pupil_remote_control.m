% pupil_remote_control.m

% (*)~----------------------------------------------------------------------------------
%  Pupil Helpers
%  Copyright (C) 2012-2016  Pupil Labs
% 
%  Distributed under the terms of the GNU Lesser General Public License (LGPL v3.0).
%  License details are in the file license.txt, distributed as part of this software.
% ----------------------------------------------------------------------------------~(*)

% Pupil Remote address
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

tic; % Measure round trip delay
zmq.core.send(socket, uint8('t'));
result = zmq.core.recv(socket,8);
fprintf('%s\n', char(result));
fprintf('Round trip command delay: %s\n', toc);

% set current Pupil time to 0.0
zmq.core.send(socket, uint8('T 0.0'));
resultZero = zmq.core.recv(socket);
fprintf('%s\n', char(resultZero));

% start recording
pause(1.0);
zmq.core.send(socket, uint8('R'));
resultR = zmq.core.recv(socket);
fprintf('Recording should start: %s\n', char(resultR));

pause(5.0);
zmq.core.send(socket, uint8('t'));
currentTime_bytes = zmq.core.recv(socket,8);
% fprintf('first measure: %s\n', char(currentTime_bytes));
currentTimeChar = char(currentTime_bytes);
currentTimeNum = str2num(currentTimeChar);
keys_start = {'topic', 'label', 'timestamp', 'duration'};
values_start = {'annotation.TrialStart', strcat(num2str(1),' started'), currentTimeNum, 1.0};
start_annotation = containers.Map(keys_start, values_start);
send_annotation(socket, start_annotation);
annotation = zmq.core.recv(socket);


pause(5.0);
zmq.core.send(socket, uint8('t'));
currentTime_bytes = zmq.core.recv(socket,8);
currentTimeChar = char(currentTime_bytes);
currentTimeNum = str2num(currentTimeChar);
keys_start = {'topic', 'label', 'timestamp', 'duration'};
values_start = {'annotation.TrialStart', strcat(num2str(1),' started'), currentTimeNum, 1.0};
start_annotation = containers.Map(keys_start, values_start);
send_annotation(socket, start_annotation);
result10 = zmq.core.recv(socket);

pause(2.0);
zmq.core.send(socket, uint8('r'));
result2 = zmq.core.recv(socket);
fprintf('Recording stopped: %s\n', char(result2));

% zmq.core.send(socket, uint8('t'));
% currentTime = zmq.core.recv(socket);


% test notification, note that you need to listen on the IPC to receive notifications!
send_notification(socket, containers.Map({'subject'}, {'calibration.should_start'}))
result3 = zmq.core.recv(socket);
fprintf('Notification received: %s\n', char(result3));

send_notification(socket, containers.Map({'subject'}, {'calibration.should_stop'}))
result3 = zmq.core.recv(socket);
fprintf('Notification received: %s\n', char(result3));


zmq.core.disconnect(socket, endpoint);
zmq.core.close(socket);

zmq.core.ctx_shutdown(ctx);
zmq.core.ctx_term(ctx);

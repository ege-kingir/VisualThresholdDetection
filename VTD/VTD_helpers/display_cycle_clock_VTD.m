function display_cycle_clock_VTD(params)
% Display clock for cycling for MAIN_GFS_Cycling.m

% Aishwarya (Ash) Bhonsle
% 22 July 2021

%% CHANGE LOG
% Following changes made by *INITIALS* made on *DATE*


%%
fprintf('\n\n\n\n\n\n Cycling warm up \n\n\n\n\n\n');

sec_deg = 0:6:354;

f = fliplr(0:1:15);
h = fliplr(16:1:59);

sec_ord = cat(2, f, h);
sec_mat = cat(1, sec_ord, sec_deg);

secondlen = 180;

Screen('TextSize', params.window, 40);
DrawFormattedText(params.window,'Please pedal at a rate of 1 lap per second.',params.centrecoords_pix(1)-400,300, [0 0 0]);

Screen('FrameOval', params.window, params.black, [params.centrecoords_pix(1)-200, params.centrecoords_pix(2)-200, params.centrecoords_pix(1)+200, params.centrecoords_pix(2)+200], 3, 3);

cyclestart = Screen('Flip', params.window);

while GetSecs - cyclestart <= params.cycledur  
    time = now;
    second = str2double(datestr(time, 'SS'));

    sec_psi = sec_mat(2, find(sec_mat(1,:) == second));        

    sec_x = secondlen * cosd(sec_psi);
    sec_y = secondlen * sind(sec_psi);

    DrawFormattedText(params.window,'Please pedal at a rate of 1 lap per second.',params.centrecoords_pix(1)-400,300, [0 0 0]);

    % Draw clock frame
    Screen('FrameOval', params.window, params.black, [params.centrecoords_pix(1)-200, params.centrecoords_pix(2)-200, params.centrecoords_pix(1)+200, params.centrecoords_pix(2)+200], 3, 3);

    Screen('DrawLine', params.window, params.black, params.centrecoords_pix(1), params.centrecoords_pix(2), params.centrecoords_pix(1) + sec_x, params.centrecoords_pix(2) - sec_y, 5); % Seconds hand

    Screen('Flip', params.window);
end;   

end
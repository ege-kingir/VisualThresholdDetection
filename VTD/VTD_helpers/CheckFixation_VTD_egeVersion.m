function ET = CheckFixation_VTD_egeVersion(ET, trial, params, output, currentCenter)
%CheckFixation CheckFixation
%   Detailed explanation goes here
if params.eyetracking == 0
    ET.IsFixation = (ET.CurrentGaze(1,1) - currentCenter(1))^2 + (ET.CurrentGaze(1,2) - currentCenter(2))^2 <= params.fixRadius^2;
elseif params.eyetracking == 1
    ET.IsFixation = ET.offset_x^2 + ET.offset_y^2 <= params.fixRadius^2;
end

if ET.IsFixation
   if ET.FixationPending
      if ET.currentTime >= ET.FixationPendingStartTime + params.TimeToEstablishFixation
         % transition from pending fixation to fixation
         ET.FixationEstablished = true;
         ET.FixationPending = false;
      end
   else
      if ET.FixationEstablished
         % Fixation and fixation is established: reset abort pending
         ET.AbortPending = false;
         
      else % new fixationPending
         ET.FixationPending = true;
         ET.FixationPendingStartTime = ET.currentTime;
      end
   end
else % no Fixation
   if ET.FixationEstablished
      if ET.AbortPending
         % Abort because of excess of blink allowance
         if ET.currentTime >= ET.AbortPendingStartTime + params.BlinkAllowance
            ET.TrialAborted = true;
            ET.AbortStartTime = ET.currentTime;
         end
      else % new AbortPending
         ET.AbortPending = true;
         ET.AbortPendingStartTime = ET.currentTime;
      end
   else
      if ET.FixationPending
         %abort a pending fixation
         ET.FixationPending = false;
      else
         if ET.currentTime >= output(trial).trialonset + params.TimeToBeginFixation
            % no fixation could be established
            ET.TrialAborted = true;
            ET.AbortStartTime = ET.currentTime;
         end
      end
   end
   
   
   
end
end


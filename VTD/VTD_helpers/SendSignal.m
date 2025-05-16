function SendSignal(obj, address, onof)
    
%      ioObj = io64;
     %address = hex2dec('EC00');%standard LPT1 output port address
%      address = hex2dec('DFB8');
     
     status = io64(obj);
     io64(obj, address, onof); 
     
     WaitSecs(0.001);
     io64(obj, address, 0); 
%      triggerout = 0;
%      triggertime = NaN; 
     
%    ioObj = io32;
%    address = hex2dec('EC00');%standard LPT1 output port address
%    status = io32(ioObj);
%    io32(ioObj, address, onof); 
end



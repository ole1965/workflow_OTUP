t = datetime('now','Format','yyyy-MM-dd''T''HHmmssSSS');S = char(t);
pointPos = findstr(gradFile,'.');
save(['rfPulses\pulseRf_',S,name,gradFile(10:pointPos-1),'.mat'],'pulseRf');

% if ~exist(['RfPulses\asInitialGuess\currentPulse_',name,'.mat'])
%     %Conversion of the complex rf pulse from rectangular to polar coordinates (by means of this the pulse is easier to constrain in the fmincon)
%     pulseSave = [abs(pulseRf(:));angle(pulseRf(:))]; 
%     iters = [];
%     zeit = [];zeit2 = datetime;
%     fcnvalues = [];
%     save(['RfPulses\asInitialGuess\currentPulse_',name,'.mat'],'iters','zeit','pulseSave','fcnvalues','zeit2');
% end
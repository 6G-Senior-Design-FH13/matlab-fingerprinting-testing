function heRangingWavGenPlot(y,cfg)
%heRangingWavGenPlot Overlays HE-LTF field names on a plot
%   heRangingWavGenPlot(Y,CFG) plots the HE Ranging waveform with overlay
%   HE-LTF symbol information for all users on first transmit
%   antenna of the first packet.
%
%   Y is the time-domain HE-LTF signal. It is a complex matrix of size
%   Ns-by-Nt where Ns represents the number of time-domain samples and Nt
%   represents the number of transmit antennas.
%
%   CFG is a format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>.

%   Copyright 2020-2022 The MathWorks, Inc.

narginchk(2,2);
if isempty(y)
    return
end

sr = wlanSampleRate(cfg.ChannelBandwidth); % Sample rate Hz
ind = double(heRangingFieldIndices(cfg,'HE-LTF'));
S = heRangingLTFInfo(cfg);

trc = wlan.internal.heTimingRelatedConstants(cfg.GuardInterval,cfg.HELTFType,4,0,0);
cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
sf = cbw*1e-3; % Scaling factor to convert bandwidth and time in ns to samples
samplesPerSym = trc.THELTFSYM*sf;

numUser = numel(cfg.User);
tick = (1/sr)*1e6; % Microseconds per sample
figure;
ax = gca;
% Plot the output of the first packet and antenna
pktInfo = cfg.validateConfig;
numFirstPktSamples = pktInfo.TxTime*cbw;
% Extract the sample required for the first packet and first antenna
y = y(1:numFirstPktSamples,1);
timeIdx = 0:tick:(length(y)-1)*tick;
plot(ax,timeIdx,20*log10(abs(y)));

hold on;
grid on;
title(sprintf('802.11az Transmission Power'))
axis([0 (length(y)-1)*tick -25 10])
xlabel('Time (us)');
ylabel('Power (dBW)');

color = 'mcrgbrkymcrgbrkymc';   
startInd = ind(1,1); % Start of HE-LTF field
% Set Y offset. This is to avoid smearing of display text
displayOffset = abs(max(y(:,1)));
userYVal = displayOffset+5;
repYVal = displayOffset+3;

% Limit display of HE-LTF field. Only display the HE-LTF field for the
% first user if the number of HE-LTF symbols is greater than 8.
if ind(numUser,2)-ind(1,1)>samplesPerSym*8
    endInd = ind(1,1)+S.NHELTFWithRepetition(1)*samplesPerSym;
    axis([0 endInd*tick -25 10]);
    numUser=1;
end

for u=1:numUser
    numRepetitions = cfg.User{u}.NumHELTFRepetition;
    plot(ax,([startInd startInd]-1)*tick,ylim(gca),'k--'); % Time index is zero-based, subtract -1 
    endInd = startInd+S.NHELTFWithoutRepetition(u)*samplesPerSym;
    updateInd = startInd:endInd;
    plot(ax,(updateInd-1)*tick,20*log10(abs(y(updateInd,1))),color(1)); hold on;

    % Display user number
    userOffset = floor(((ind(u,2)+1-ind(u,1))/2+ind(u,1))/sr*1e6); % Offset in us
    text(ax,userOffset,userYVal,sprintf('User-%d',u),'FontSize',8,'HorizontalAlignment','center')

    % Display HE-LTF repetition
    repLength = (ind(u,2)+1-ind(u,1))/numRepetitions;
    repOffset = floor((repLength/2+ind(u,1)+(0:repLength:(numRepetitions-1)*repLength))/sr*1e6); % Offset in us

    % First repetition
    text(ax,repOffset(1),repYVal,sprintf('Rep-%d',1),'FontSize',8,'HorizontalAlignment','center')
    % Remaining repetitions
    for r=2:cfg.User{u}.NumHELTFRepetition
        text(ax,repOffset(r),repYVal,sprintf('Rep-%d',r),'FontSize',8,'HorizontalAlignment','center')
        yAxis = ax.YLim(1)/2:ax.YLim(2)/2;
        plot(ax,(endInd-1)*tick*ones(numel(yAxis),1),yAxis,'k--'); hold on;
        startInd = endInd;
        endInd = startInd+S.NHELTFWithoutRepetition(u)*samplesPerSym;
        updateInd = startInd:endInd;
        plot(ax,(updateInd-1)*tick,20*log10(abs(y(updateInd,1))),color(r)); hold on;
    end 

    % Plot the last repetition of HE-LTF symbol
    plot(ax,([endInd endInd]-1)*tick,ylim(gca),'k--');
    startInd = endInd;
end

end

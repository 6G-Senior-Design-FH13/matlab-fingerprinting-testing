%Jean Lee's Code
TR = stlread("office.stl");
scale = 0.9;
scaledPts = TR.Points * scale;
TR_scaled = triangulation(TR.ConnectivityList,scaledPts);
viewer = siteviewer('SceneModel',TR_scaled);

% Define transmitter
tx = txsite("cartesian", ...
    "AntennaPosition",[0.01; 1; 2.5], ...
    "TransmitterFrequency",1e11); %100GHz max
show(tx,"ShowAntennaHeight",false) 

% Define receiver
rx = rxsite("cartesian", ...
    "AntennaPosition",[4.5; 3.5; 0.8]);
show(rx,"ShowAntennaHeight",false)

txArraySize = [1 1];
rxArraySize = [1 1];

los(tx,rx); %line of sight

pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","image", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","wood"); 

% Perform ray tracing
ray = raytrace(tx,rx,pm);
% ray = ray{1};
plot(ray{1});

snrs = [10 15 20];  % range of SNR --> different noise conditions
chanBW = "CBW40"; 
cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
    "NumTransmitAntennas",prod(txArraySize), ...
    "SecureHELTF",false, ...
    "GuardInterval",1.6);

[features,labels] = dlPositioningGenerateDataSet(ray,rx,tx,cfg,snrs,txArraySize,rxArraySize);

function [cir, labels] = dlPositioningGenerateDataSet(rays, STAs, APs, cfg, snrs, txArraySize, rxArraySize)
%HELPERGENERATEDATASET Generate a Dataset of 802.11az CIR Fingerprints 
%   dlPositioningGenerateDataSet(RAYS, STAS, APS, CFG, SNRS) create
%   multiple channels from multi-path propagation objects, RAYS, between
%   all transmitters, APs, and receivers, STAs. Generates 802.11az waveform
%   to be passed parameterized by, CFG, and sets the range of awgn noise
%   values, SNRS, to be added to the received signal. Returns a matrix with
%   the CIR of all channel realizations and their positions and locations,
%   LABELS.

%   Copyright 2020-2021 The MathWorks, Inc.

ofdmSymbolOffset = 0.75;

numChan = numel(rays);
txWaveform = single(heRangingWaveformGenerator(cfg)); % Generate 802.11az packet

ofdmInfo = wlanHEOFDMInfo('HE-LTF',cfg.ChannelBandwidth,cfg.GuardInterval);

% Create empty cir (feature) and loc (labels) matrices. cir is of size Ns x
% Nsts*Nr x Nsnr x Ntx*Nrx 4-D matrix This is reshaped to Ns x Nsts*Nr x
% Naps x Nsnr*Nstas for easier processing and use in the CNN.
cir = zeros([ofdmSymbolOffset*ofdmInfo.CPLength prod(txArraySize)*prod(rxArraySize) length(snrs) numChan],'single');
labels.position = zeros([3 numChan]);
for i = 1:numChan
    txn = mod(i-1,height(rays))+1;
    rxn = ceil(i/height(rays));
    if isempty(rays{i})
        % If no rays were received from a tx at this rx 0 the
        % indices of the matrix for that location and store the
        % position.
        cir(:,:,:,i) = 0;    
    else
        % Generates the channel estimate/returns the CIR.
        cir(:,:,:,i) = generateCIR(rays{i},APs(txn),STAs(rxn),cfg,txWaveform,ofdmInfo,snrs,ofdmSymbolOffset,txArraySize,rxArraySize);     
    end       
    labels.position(:,i) = [STAs(rxn).AntennaPosition];
    labels.class(i) = categorical(cellstr(STAs(rxn).Name));  

    % Displays progress (10% intervals)
    if mod(i,floor(numChan/10))==0
        qt = ceil(i/(numChan/10));
        disp(['Generating Dataset: ', num2str(10*qt), '% complete.'])
    end
end

% Ns x Nsts*Nr x Nsnr x Naps x Nstas
cir = reshape(cir,[ofdmSymbolOffset*ofdmInfo.CPLength prod(txArraySize)*prod(rxArraySize) length(snrs) numel(APs) numel(STAs)]);
% Ns x Nsts*Nr x Naps x Nsnr x Nstas
cir = permute(cir,[1 2 4 3 5]);
% Ns x Nsts*Nr x Naps x Nsnr*Nstas
cir = reshape(cir,[size(cir,1) size(cir,2) size(cir,3) size(cir,4)*size(cir,5)]);

labels.position = labels.position(:, 1:height(rays):end);
labels.class = labels.class(:, 1:height(rays):end); % Remove duplicated locations

% Create and scale training labels from rx locations to correct size
labels.position = repelem(labels.position, 1, length(snrs));
labels.class = repelem(labels.class, 1, length(snrs));

end

function cir = generateCIR(rays,AP,STA,cfg,tx,ofdmInfo,snr,ofdmSymbolOffset,txArraySize,rxArraySize)
%GENERATECIR Generate CIR fingerprint for MIMO channel.
%   generateCIR(RAYS,AP,STA,CFG,TX,OFDMINFO,SNR) returns a CIR fingerprint
%   by constructing a channel from multi-path propagation objects, RAYS,
%   between a single AP and STA. The configuration, CFG, and OFDMINFO are
%   used to set channel parameters. TX is a packet to be passed through the
%   channel and SNR defines the noise to be added at the receiver.
 
    rtChan = comm.RayTracingChannel(rays,AP,STA); % Create channel
    rtChan.SampleRate = wlanSampleRate(cfg.ChannelBandwidth);
    rtChan.ReceiverVirtualVelocity = [0; 0; 0]; % Stationary Receiver
    rtChan.NormalizeChannelOutputs = false;
    
    rxChan = rtChan(tx); % Pass waveform through channel
    
    % Create matrix for CIR to be stored in. 
    % Dimensions are Ns x Nsts*Nr x Nsnr
    cir = zeros([ofdmSymbolOffset*ofdmInfo.CPLength prod(txArraySize)*prod(rxArraySize) length(snr)],'single');
    
    % Adjust power of noise added such that the SNR is per active
    % subcarrier
    snrAdj = snr-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
    
    % Add noise to the the received waveform for each snr value,
    % Perform synchronization, channel estimation and extract the CIR for
    % each.
    for i=1:length(snr)
        rx = awgn(rxChan,snrAdj(i)); 
        chanEst = heRangingSynchronize(double(rx),cfg); % Perform synchronization and channel estimation.
        if isempty(chanEst)
            continue % Synchronization fails
        end
        
        % Trim CIR to make data more manageable. Assume the CIR fits into
        % the useful portion of the CP (otherwise ISI present)
        cirRaw = helperChannelImpulseResponse(single(chanEst),ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices);
        cir(:,:,i) = reshape(abs(cirRaw(1:ofdmSymbolOffset*ofdmInfo.CPLength,:,:)),ofdmSymbolOffset*ofdmInfo.CPLength,[]);
    end
end








% 5G waveform.....


% carrier = nrSCSCarrierConfig('NSizeGrid',100);
% bwp = nrWavegenBWPConfig('NStartBWP',carrier.NStartGrid+10);
% ssb = nrWavegenSSBurstConfig('BlockPattern','Case A');
% pdcch = nrWavegenPDCCHConfig('AggregationLevel',2,'AllocatedCandidate',4);
% coreset = nrCORESETConfig;
% coreset.FrequencyResources = [1 1 1 1];
% coreset.Duration = 3;
% ss = nrSearchSpaceConfig;
% ss.NumCandidates = [8 4 0 0 0];
% pdsch = nrWavegenPDSCHConfig( ...
%     'Modulation','16QAM','TargetCodeRate',658/1024,'EnablePTRS',true);
% dmrs = nrPDSCHDMRSConfig('DMRSTypeAPosition',3);
% pdsch.DMRS = dmrs;
% ptrs = nrPDSCHPTRSConfig('TimeDensity',2);
% pdsch.PTRS = ptrs;
% csirs = nrWavegenCSIRSConfig('RowNumber',4,'RBOffset',10);
% cfgDL = nrDLCarrierConfig( ...
%     'FrequencyRange','FR1', ...
%     'ChannelBandwidth',40, ...
%     'NumSubframes',20, ...
%     'SCSCarriers',{carrier}, ...
%     'BandwidthParts',{bwp}, ...
%     'SSBurst',ssb, ...
%     'CORESET',{coreset}, ...
%     'SearchSpaces',{ss}, ...
%     'PDCCH',{pdcch}, ...
%     'PDSCH',{pdsch}, ...
%     'CSIRS',{csirs});
% 
% snrs = [10 15 20]; 
% 
% txwaveform = nrWaveformGenerator(cfgDL); %Generate 5G downlink for a single user
% 
% numChan = numel(r);
% 
% ofdmSymbolOffset = 0.75;
% 
% ofdmInfo = nrOFDMInfo(carrier);
% 
% cir = zeros([ofdmSymbolOffset*ofdmInfo.CyclicPrefixLengths prod(txArraySize)*prod(rxArraySize) length(snrs) numChan],'single');
% labels.position = zeros([3 numChan]);
% 
% 
% 
% % 
% % SNR = 1;
% % N0 = 1/sqrt(2.0*rxArraySize*double(odfmInfo.Nfft)*SNR);
% % w = (1/sqrt(pdsch.NumLayers))*ones(pdsch.NumLayers,txArraySize);
% % [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
% % pdschBits = randi([0 1],pdschInfo.G,1);
% % pdschSymbols = nrPDSCH(carrier,pdsch,pdschBits);
% % pdschSymbolsPrecoded = pdschSymbols*w;
% % 
% % pdschGrid = nrResourceGrid(carrier,txArraySize);
% % [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
% % pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;
% % txWaveform = nrOFDMModulate(carrier,pdschGrid);
% % rxWaveform = txWaveform/sqrt(rxArraySize);
% % rxNoise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
% % % OFDM demodulation
% % rxSignalGrid = nrOFDMDemodulate(carrier,rxWaveform);
% % rxNoiseGrid = nrOFDMDemodulate(carrier,rxNoise);
% % 
% % % PDSCH symbols extraction
% % rxPDSCHSymbols = rxSignalGrid(pdschAntIndices);
% % 
% % Sre = (1/waveformInfo.Nfft.^2)*rms(rxPDSCHSymbols).^2;
% % Nre = (1/waveformInfo.Nfft)*rms(rxNoise).^2;
% 
% 
% 
% for i = 1:numChan
%     txn = mod(i-1,height(r))+1;
%     rxn = ceil(i/height(r));
%     if isempty(r{i})
%         % If no rays were received from a tx at this rx 0 the
%         % indices of the matrix for that location and store the
%         % position.
%         cir(:,:,:,i) = 0;    
%     else
%         % Generates the channel estimate/returns the CIR.
%         cir(:,:,:,i) = generateCIR(r{i},tx(txn),rx(rxn),cfgDL,txWaveform,ofdmInfo,snrs,ofdmSymbolOffset, txArraySize, rxArraySize);     
%     end       
%     labels.position(:,i) = [rx(rxn).AntennaPosition];
%     labels.class(i) = categorical(cellstr(rx(rxn).Name));  
% 
%     % Displays progress (10% intervals)
%     if mod(i,floor(numChan/10))==0
%         qt = ceil(i/(numChan/10));
%         disp(['Generating Dataset: ', num2str(10*qt), '% complete.'])
%     end
% end
% 
% 
% 
% % Ns x Nsts*Nr x Nsnr x Naps x Nstas
% cir = reshape(cir,[ofdmSymbolOffset*ofdmInfo.CyclicPrefixLengths prod(txArraySize)*prod(rxArraySize) length(snrs) numel(tx) numel(rx)]);
% % Ns x Nsts*Nr x Naps x Nsnr x Nstas
% cir = permute(cir,[1 2 4 3 5]);
% % Ns x Nsts*Nr x Naps x Nsnr*Nstas
% cir = reshape(cir,[size(cir,1) size(cir,2) size(cir,3) size(cir,4)*size(cir,5)]);
% 
% labels.position = labels.position(:, 1:height(r):end);
% labels.class = labels.class(:, 1:height(r):end); % Remove duplicated locations
% 
% % Create and scale training labels from rx locations to correct size
% labels.position = repelem(labels.position, 1, length(snrs));
% labels.class = repelem(labels.class, 1, length(snrs));
% 
% 
% 

% function cir = generateCIR(rays,AP,STA,cfg,waveform,ofdmInfo,snr,ofdmSymbolOffset, txArraySize, rxArraySize)
% %GENERATECIR Generate CIR fingerprint for MIMO channel.
% %   generateCIR(RAYS,AP,STA,CFG,TX,OFDMINFO,SNR) returns a CIR fingerprint
% %   by constructing a channel from multi-path propagation objects, RAYS,
% %   between a single AP and STA. The configuration, CFG, and OFDMINFO are
% %   used to set channel parameters. TX is a packet to be passed through the
% %   channel and SNR defines the noise to be added at the receiver.
%  
%     rtChan = comm.RayTracingChannel(rays,AP,STA); % Create channel
%     rtChan.SampleRate = wlanSampleRate(cfg.ChannelBandwidth);
%     rtChan.ReceiverVirtualVelocity = [0; 0; 0]; % Stationary Receiver
%     rtChan.NormalizeChannelOutputs = false;
%     
%     rxChan = rtChan(waveform); % Pass waveform through channel
%     
%     % Create matrix for CIR to be stored in. 
%     % Dimensions are Ns x Nsts*Nr x Nsnr
%     cir = zeros([ofdmSymbolOffset*ofdmInfo.CyclicPrefixLengths prod(txArraySize)*prod(rxArraySize) length(snr)],'single');
%     
%     % Adjust power of noise added such that the SNR is per active
%     % subcarrier
% 
%     channel = nrTDLChannel;
%     channel.DelayProfile = "TDL-C";
%     channel.NumTransmitAntennas = txArraySize;
%     channel.NumReceiveAntennas = rxArraySize;
%     
%     [txWaveform,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);
%     chInfo = info(channel);
%     maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
%     txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];
%     
%     [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
%     
% 
%     
%     % noise = generateAWGN(SNRdB,nRxAnts,waveformInfo.Nfft,size(rxWaveform));
%     rxWaveform = rxWaveform + noise;
%     
%     if perfectEstimation
%         % Get path filters for perfect timing estimation
%         pathFilters = getPathFilters(channel); 
%         [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
%     else
%         [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
%         offset = hSkipWeakTimingOffset(offset,t,mag);
%     end
%     rxWaveform = rxWaveform(1+offset:end,:);
%     
% 
% 
% 
% 
% 
% 
%     snrAdj = snr-pow2db(Sre/Nre);
% 
%     % Add noise to the the received waveform for each snr value,
%     % Perform synchronization, channel estimation and extract the CIR for
%     % each.
%     for i=1:length(snr)
%         SNR = 10^(SNRdB/10); % Calculate linear noise gain
%         N0 = 1/sqrt(2.0*nRxAnts*double(Nfft)*SNR);
%         noise = N0*complex(randn(sizeRxWaveform),randn(sizeRxWaveform));
%         rxWavefore = rxWaveform + noise;
% %         rx = awgn(rxChan,snrAdj(i)); 
%         chanEst = heRangingSynchronize(double(rx),cfg); % Perform synchronization and channel estimation.
%         if isempty(chanEst)
%             continue % Synchronization fails
%         end
%         
%         % Trim CIR to make data more manageable. Assume the CIR fits into
%         % the useful portion of the CP (otherwise ISI present)
%         cirRaw = helperChannelImpulseResponse(single(chanEst),ofdmInfo.NFFT,ofdmInfo.CyclicPrefixLengths,ofdmInfo.ActiveFFTIndices);
%         cir(:,:,i) = reshape(abs(cirRaw(1:ofdmSymbolOffset*ofdmInfo.CyclicPrefixLengths,:,:)),ofdmSymbolOffset*ofdmInfo.CyclicPrefixLengths,[]);
%     end
% end


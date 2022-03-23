% % function [cir, labels] = dlPositioningGenerateDataSet(rays, STAs, APs, cfg, snrs)
% % %HELPERGENERATEDATASET Generate a Dataset of 802.11az CIR Fingerprints 
% % %   dlPositioningGenerateDataSet(RAYS, STAS, APS, CFG, SNRS) create
% % %   multiple channels from multi-path propagation objects, RAYS, between
% % %   all transmitters, APs, and receivers, STAs. Generates 802.11az waveform
% % %   to be passed parameterized by, CFG, and sets the range of awgn noise
% % %   values, SNRS, to be added to the received signal. Returns a matrix with
% % %   the CIR of all channel realizations and their positions and locations,
% % %   LABELS.
% % 
% % %   Copyright 2020-2021 The MathWorks, Inc.
% % 
% % ofdmSymbolOffset = 0.75;
% % 
% % numChan = numel(rays);
% % txWaveform = single(heRangingWaveformGenerator(cfg)); % Generate 802.11az packet
% % 
% % ofdmInfo = wlanHEOFDMInfo('HE-LTF',cfg.ChannelBandwidth,cfg.GuardInterval);
% % 
% % % Create empty cir (feature) and loc (labels) matrices. cir is of size Ns x
% % % Nsts*Nr x Nsnr x Ntx*Nrx 4-D matrix This is reshaped to Ns x Nsts*Nr x
% % % Naps x Nsnr*Nstas for easier processing and use in the CNN.
% % cir = zeros([ofdmSymbolOffset*ofdmInfo.CPLength cfg.User{1}.NumSpaceTimeStreams*prod(STAs(1).Antenna.Size) length(snrs) numChan],'single');
% % labels.position = zeros([3 numChan]);
% % for i = 1:numChan
% %     txn = mod(i-1,height(rays))+1;
% %     rxn = ceil(i/height(rays));
% %     if isempty(rays{i})
% %         % If no rays were received from a tx at this rx 0 the
% %         % indices of the matrix for that location and store the
% %         % position.
% %         cir(:,:,:,i) = 0;    
% %     else
% %         % Generates the channel estimate/returns the CIR.
% %         cir(:,:,:,i) = generateCIR(rays{i},APs(txn),STAs(rxn),cfg,txWaveform,ofdmInfo,snrs,ofdmSymbolOffset);     
% %     end       
% %     labels.position(:,i) = [STAs(rxn).AntennaPosition];
% %     labels.class(i) = categorical(cellstr(STAs(rxn).Name));  
% % 
% %     % Displays progress (10% intervals)
% %     if mod(i,floor(numChan/10))==0
% %         qt = ceil(i/(numChan/10));
% %         disp(['Generating Dataset: ', num2str(10*qt), '% complete.'])
% %     end
% % end
% % 
% % % Ns x Nsts*Nr x Nsnr x Naps x Nstas
% % cir = reshape(cir,[ofdmSymbolOffset*ofdmInfo.CPLength cfg.User{1}.NumSpaceTimeStreams*prod(STAs(1).Antenna.Size) length(snrs) numel(APs) numel(STAs)]);
% % % Ns x Nsts*Nr x Naps x Nsnr x Nstas
% % cir = permute(cir,[1 2 4 3 5]);
% % % Ns x Nsts*Nr x Naps x Nsnr*Nstas
% % cir = reshape(cir,[size(cir,1) size(cir,2) size(cir,3) size(cir,4)*size(cir,5)]);
% % 
% % labels.position = labels.position(:, 1:height(rays):end);
% % labels.class = labels.class(:, 1:height(rays):end); % Remove duplicated locations
% % 
% % % Create and scale training labels from rx locations to correct size
% % labels.position = repelem(labels.position, 1, length(snrs));
% % labels.class = repelem(labels.class, 1, length(snrs));
% % 
% % end
% % 
% % function cir = generateCIR(rays,AP,STA,cfg,tx,ofdmInfo,snr,ofdmSymbolOffset)
% % %GENERATECIR Generate CIR fingerprint for MIMO channel.
% % %   generateCIR(RAYS,AP,STA,CFG,TX,OFDMINFO,SNR) returns a CIR fingerprint
% % %   by constructing a channel from multi-path propagation objects, RAYS,
% % %   between a single AP and STA. The configuration, CFG, and OFDMINFO are
% % %   used to set channel parameters. TX is a packet to be passed through the
% % %   channel and SNR defines the noise to be added at the receiver.
% %  
% %     rtChan = comm.RayTracingChannel(rays,AP,STA); % Create channel
% %     rtChan.SampleRate = wlanSampleRate(cfg.ChannelBandwidth);
% %     rtChan.ReceiverVirtualVelocity = [0; 0; 0]; % Stationary Receiver
% %     rtChan.NormalizeChannelOutputs = false;
% %     
% %     rxChan = rtChan(tx); % Pass waveform through channel
% %     
% %     % Create matrix for CIR to be stored in. 
% %     % Dimensions are Ns x Nsts*Nr x Nsnr
% %     cir = zeros([ofdmSymbolOffset*ofdmInfo.CPLength cfg.User{1}.NumSpaceTimeStreams*prod(STA.Antenna.Size) length(snr)],'single');
% %     
% %     % Adjust power of noise added such that the SNR is per active
% %     % subcarrier
% %     snrAdj = snr-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
% %     
% %     % Add noise to the the received waveform for each snr value,
% %     % Perform synchronization, channel estimation and extract the CIR for
% %     % each.
% %     for i=1:length(snr)
% %         rx = awgn(rxChan,snrAdj(i)); 
% %         chanEst = heRangingSynchronize(double(rx),cfg); % Perform synchronization and channel estimation.
% %         if isempty(chanEst)
% %             continue % Synchronization fails
% %         end
% %         
% %         % Trim CIR to make data more manageable. Assume the CIR fits into
% %         % the useful portion of the CP (otherwise ISI present)
% %         cirRaw = helperChannelImpulseResponse(single(chanEst),ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices);
% %         cir(:,:,i) = reshape(abs(cirRaw(1:ofdmSymbolOffset*ofdmInfo.CPLength,:,:)),ofdmSymbolOffset*ofdmInfo.CPLength,[]);
% %     end
% % end
% 
% function [cir] = dlPositioningGenerateDataSet(rays, STAs, APs, cfg, snrs)
% %HELPERGENERATEDATASET Generate a Dataset of 802.11az CIR Fingerprints 
% %   dlPositioningGenerateDataSet(RAYS, STAS, APS, CFG, SNRS) create
% %   multiple channels from multi-path propagation objects, RAYS, between
% %   all transmitters, APs, and receivers, STAs. Generates 802.11az waveform
% %   to be passed parameterized by, CFG, and sets the range of awgn noise
% %   values, SNRS, to be added to the received signal. Returns a matrix with
% %   the CIR of all channel realizations and their positions and locations,
% %   LABELS.
% 
% %   Copyright 2020-2021 The MathWorks, Inc.
% 
% ofdmSymbolOffset = 0.75;
% 
% numChan = numel(rays);
% txWaveform = single(heRangingWaveformGenerator(cfg)); % Generate 802.11az packet
% 
% ofdmInfo = wlanHEOFDMInfo('HE-LTF',cfg.ChannelBandwidth,cfg.GuardInterval);
% 
% % Create empty cir (feature) and loc (labels) matrices. cir is of size Ns x
% % Nsts*Nr x Nsnr x Ntx*Nrx 4-D matrix This is reshaped to Ns x Nsts*Nr x
% % Naps x Nsnr*Nstas for easier processing and use in the CNN.
% cir = zeros([ofdmSymbolOffset*ofdmInfo.CPLength cfg.User{1}.NumSpaceTimeStreams length(snrs) numChan],'single');
% %labels.position = zeros([3 numChan]);
% for i = 1:numChan
%     txn = mod(i-1,height(rays))+1;
%     rxn = ceil(i/height(rays));
%     if isempty(rays{i})
%         % If no rays were received from a tx at this rx 0 the
%         % indices of the matrix for that location and store the
%         % position.
%         cir(:,:,:,i) = 0;    
%     else
%         % Generates the channel estimate/returns the CIR.
%         cir(:,:,:,i) = generateCIR(rays{i},APs(txn),STAs(rxn),cfg,txWaveform,ofdmInfo,snrs,ofdmSymbolOffset);     
%     end       
%     %labels.position(:,i) = [STAs(rxn).AntennaPosition];
%     %labels.class(i) = categorical(cellstr(STAs(rxn).Name));  
% 
%     % Displays progress (10% intervals)
%     if mod(i,floor(numChan/10))==0
%         qt = ceil(i/(numChan/10));
%         disp(['Generating Dataset: ', num2str(10*qt), '% complete.'])
%     end
% end
% 
% % Ns x Nsts*Nr x Nsnr x Naps x Nstas
% cir = reshape(cir,[ofdmSymbolOffset*ofdmInfo.CPLength cfg.User{1}.NumSpaceTimeStreams length(snrs) numel(APs) numel(STAs)]);
% % Ns x Nsts*Nr x Naps x Nsnr x Nstas
% cir = permute(cir,[1 2 4 3 5]);
% % Ns x Nsts*Nr x Naps x Nsnr*Nstas
% cir = reshape(cir,[size(cir,1) size(cir,2) size(cir,3) size(cir,4)*size(cir,5)]);
% 
% %labels.position = labels.position(:, 1:height(rays):end);
% %labels.class = labels.class(:, 1:height(rays):end); % Remove duplicated locations
% 
% % Create and scale training labels from rx locations to correct size
% %labels.position = repelem(labels.position, 1, length(snrs));
% %labels.class = repelem(labels.class, 1, length(snrs));
% 
% end
% 
% function cir = generateCIR(rays,AP,STA,cfg,tx,ofdmInfo,snr,ofdmSymbolOffset)
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
%     rxChan = rtChan(tx); % Pass waveform through channel
%     
%     % Create matrix for CIR to be stored in. 
%     % Dimensions are Ns x Nsts*Nr x Nsnr
%     cir = zeros([ofdmSymbolOffset*ofdmInfo.CPLength cfg.User{1}.NumSpaceTimeStreams length(snr)],'single');
%     
%     % Adjust power of noise added such that the SNR is per active
%     % subcarrier
%     snrAdj = snr-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
%     
%     % Add noise to the the received waveform for each snr value,
%     % Perform synchronization, channel estimation and extract the CIR for
%     % each.
%     for i=1:length(snr)
%         rx = awgn(rxChan,snrAdj(i)); 
%         chanEst = heRangingSynchronize(double(rx),cfg); % Perform synchronization and channel estimation.
%         if isempty(chanEst)
%             continue % Synchronization fails
%         end
%         
%         % Trim CIR to make data more manageable. Assume the CIR fits into
%         % the useful portion of the CP (otherwise ISI present)
%         cirRaw = helperChannelImpulseResponse(single(chanEst),ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices);
%         cir(:,:,i) = reshape(abs(cirRaw(1:ofdmSymbolOffset*ofdmInfo.CPLength,:,:)),ofdmSymbolOffset*ofdmInfo.CPLength,[]);
%     end
% end



function [range, chanEst] = dlPositioningGenerateDataSet(rays, STAs, APs, cfg, snrs)
    % Adjust power of noise added such that the SNR is per active
    % subcarrier
%     snrAdj = snr-10*log10();
    numChan = numel(rays);

% % % % % % % % % % % % % % % % %     
    sr = 40e6;
    t = 0:1/sr:159*(1/sr);
    peField = chirp(t,1e6,t(end),20e6)'; % 20e6 = sr/2

    
    for i = 1:numChan
        txn = mod(i-1,height(rays))+1;
        rxn = ceil(i/height(rays));


        txWaveform = single(heRangingWaveformGenerator(cfg)); % Generate 802.11az packet
%         figure(1);
%         pspectrum(txWaveform,40e6,'spectrogram', 'TimeResolution',0.000001, 'OverlapPercent',99);


        rtChan = comm.RayTracingChannel(rays{i},APs(txn),STAs(rxn)); % Create channel
        rtChan.SampleRate = wlanSampleRate(cfg.ChannelBandwidth);
        rtChan.ReceiverVirtualVelocity = [0; 0; 0]; % Stationary Receiver
        rtChan.NormalizeChannelOutputs = false;
        
        rxChan = rtChan(txWaveform); % Pass waveform through channel
        
        % Add noise to the the received waveform for each snr value,
        % Perform synchronization, channel estimation
        for j=1:length(snrs)
            rx = awgn(rxChan,snrs(j)); 
            chanEst = heRangingSynchronize(double(rx),cfg); % Perform synchronization and channel estimation.
            if isempty(chanEst)
                continue % Synchronization fails
            end
        end
    end 
        r = xcorr(chanEst, txWaveform);
        L = 160;
        Fs = sr;
        Y = fft(r);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        figure(11);
        plot(f,P1);
        hold on;
        bw = 40e6; %CBW40 = 40MHz
        T = 3.5e-6;
        [maxY, indexMaxY] = max(P1);
        fp = f(indexMaxY(1));
        range = bw/(T*fp);
%         disp(range);
end

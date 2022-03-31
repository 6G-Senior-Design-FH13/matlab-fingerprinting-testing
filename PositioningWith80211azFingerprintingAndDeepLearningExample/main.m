% mapFileName = "office.stl";
% viewer = siteviewer("SceneModel",mapFileName,"Transparency",0.25);
% S = RandStream("mt19937ar","Seed",5489); % Set the RNG for reproducibility.
% RandStream.setGlobalStream(S);
% 
% txArraySize = [4 1]; % Linear transmit array
% rxArraySize = [1 1]; % Linear receive array
% chanBW = "CBW40"; 
% distribution = "uniform";
% staSeparation = 1; % STA separation, in meters, used only when the distribution is uniform
% numSTAs = 50;       % Number of STAs, used only when distribution is random
% task = "localization"; 
% if distribution == "uniform"
%     [APs,STAs] = dlPositioningCreateEnvironment(txArraySize,rxArraySize,staSeparation,"uniform");
% else
%     [APs,STAs] = dlPositioningCreateEnvironment(txArraySize,rxArraySize,numSTAs,"random");
% end
% show(APs)
% show(STAs,'ShowAntennaHeight',false,'IconSize',[16 16]);
% pm = propagationModel("raytracing", ...
%     "CoordinateSystem","cartesian", ...
%     "SurfaceMaterial","wood", ...
%     "MaxNumReflections",1);
% % Perform ray tracing for all transmitters and receivers in parallel
% rays = raytrace(APs,STAs,pm,"Map",mapFileName);
% size(rays)
% hide(STAs);
% show(STAs(30),'IconSize',[32 32]);
% plot([rays{:,30}],'ColorLimits',[50 95]);
% snrs = [10 15 20]; 
% cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
%     "NumTransmitAntennas",prod(txArraySize), ...
%     "SecureHELTF",false, ...
%     "GuardInterval",1.6);
% cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);
% [features,labels] = dlPositioningGenerateDataSet(rays,STAs,APs,cfg,snrs);
close all; clearvars; clc;




mapFileName = "office.stl";
% viewer = siteviewer("SceneModel",mapFileName,"Transparency",0.25);
distribution = "uniform";
txArraySize = [1 1]; % Linear transmit array
rxArraySize = [1 1]; % Linear receive array
chanBW = "CBW40"; 
staSeparation = .5; % STA separation, in meters, used only when the distribution is uniform
numSTAs = 480;  
S = RandStream("mt19937ar","Seed",5489); % Set the RNG for reproducibility.
RandStream.setGlobalStream(S);
if distribution == "uniform"
    [APs,STAs] = dlPositioningCreateEnvironment(txArraySize,rxArraySize,staSeparation,"uniform");
else
    [APs,STAs] = dlPositioningCreateEnvironment(txArraySize,rxArraySize,numSTAs,"random");
end
%show(APs)
%show(STAs,'ShowAntennaHeight',false,'IconSize',[16 16]);
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "SurfaceMaterial","wood", ...
    "MaxNumReflections",3);
rays = raytrace(APs,STAs,pm,"Map",mapFileName);
snr = 10; 
cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
    "NumTransmitAntennas",prod(txArraySize), ...
    "SecureHELTF",false, ...
    "GuardInterval",1.6);
cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);
[p] = dlPositioningGenerateDataSet(rays,STAs,APs,cfg,snr);

% 
% mapFileName = "office.stl";
% TR = stlread("office.stl");
% scale = 0.9;
% scaledPts = TR.Points * scale;
% TR_scaled = triangulation(TR.ConnectivityList,scaledPts);
% 
% viewer = siteviewer('SceneModel',TR_scaled);
% 
% xrx = [0.1 4.4];
% yrx = [0.1 6.9];
% zrx = [2.1];
% antPosAP = [kron(xrx, ones(1, length(yrx))); ...
%           repmat(yrx, 1, length(xrx)); ...
%           zrx*ones(1, length(xrx)*length(yrx))];
% % antPosAP = [0.1 0.1 4.4 4.4;
% %             0.1 6.9 0.1 6.9;
% %             1.1 2.1 1.1 2.1];
% 
% snrs = 10; 
% chanBW = "CBW40"; 
% txArraySize = [1 1];
% rxArraySize = [1 1];
% % % % % % % % % % % % % % % % % % 
% fc = 1e9; % Set the carrier frequency (Hz) (support for ray tracing -- 100MHz to 1GHz)
% % % % % % % % % % % % % % % % % % % 
% lambda = physconst("lightspeed")/fc;
% 
% txArray = arrayConfig("Size", txArraySize, "ElementSpacing", 2*lambda);
% 
% tx = txsite("cartesian", ...
%     "Antenna", txArray, ...
%     "AntennaPosition", antPosAP,...
%     "TransmitterFrequency", fc);
% show(tx,"ShowAntennaHeight",false);
% rx = rxsite("cartesian", ...
%     "AntennaPosition",[1; 2; 2]);
% show(rx,"ShowAntennaHeight",false)
% 
% pm = propagationModel("raytracing", ...
%     "CoordinateSystem","cartesian", ...
%     "Method","image", ...
%     "MaxNumReflections",2, ...
%     "SurfaceMaterial","wood"); 
% % r = raytrace(tx,rx,pm);
% % plot(r{1});
% 
% rays1 = raytrace(tx,rx,pm,"Map",mapFileName);
% 
% cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
%     "NumTransmitAntennas",1, ...
%     "SecureHELTF",false, ...
%     "GuardInterval",1.6);
% cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);
% 
% % y = dlPositioningGenerateDataSet(rays,rx,tx,cfg,snrs);
% [r1] = dlPositioningGenerateDataSet(rays1,rx,tx,cfg,snrs);
% 
% % figure(3);
% % pspectrum(y1,40e6,'spectrogram', 'TimeResolution',0.000001, 'OverlapPercent',99);
% 
% 
% 

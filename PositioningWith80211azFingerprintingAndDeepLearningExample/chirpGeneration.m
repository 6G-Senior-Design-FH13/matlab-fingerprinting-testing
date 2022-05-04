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


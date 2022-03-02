mapFileName = "office.stl";
viewer = siteviewer("SceneModel",mapFileName,"Transparency",0.25);
distribution = "uniform";
txArraySize = [4 1]; % Linear transmit array
rxArraySize = [1 1]; % Linear transmit array
chanBW = "CBW40"; 
staSeparation = .5; % STA separation, in meters, used only when the distribution is uniform
numSTAs = 300;  
S = RandStream("mt19937ar","Seed",5489); % Set the RNG for reproducibility.
RandStream.setGlobalStream(S);
if distribution == "uniform"
    [APs,STAs] = DEFdlPositioningCreateEnvironment(txArraySize,rxArraySize,staSeparation,"uniform");
else
    [APs,STAs] = DEFdlPositioningCreateEnvironment(txArraySize,rxArraySize,numSTAs,"random");
end
show(APs)
show(STAs,'ShowAntennaHeight',false,'IconSize',[16 16]);
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "SurfaceMaterial","wood", ...
    "MaxNumReflections",1);
rays = raytrace(APs,STAs,pm,"Map",mapFileName);
snr = 10; 
cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
    "NumTransmitAntennas",prod(txArraySize), ...
    "SecureHELTF",false, ...
    "GuardInterval",1.6);
cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);
[features,labels] = DEFdlPositioningGenerateDataSet(rays,STAs,APs,cfg,snr);

if distribution == "uniform"
    save('output/feats4T.mat', 'features', '-mat' );
    save('output/labels4T.mat', 'lp', '-mat' );
else
    save('output/feats4TValidate.mat', 'features', '-mat' );
    save('output/labels4TValidate.mat', 'lp', '-mat' );
end

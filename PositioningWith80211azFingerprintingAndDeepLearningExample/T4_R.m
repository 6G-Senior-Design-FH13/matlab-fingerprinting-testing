mapFileName = "office.stl";
viewer = siteviewer("SceneModel",mapFileName,"Transparency",0.25);
txArraySize = [4 1]; % Linear transmit array
rxArraySize = [4 1]; % Linear transmit array
chanBW = "CBW40"; 
staSeparation = .5; % STA separation, in meters, used only when the distribution is uniform
[APs,STAs] = DEFdlPositioningCreateEnvironment(txArraySize,rxArraySize,staSeparation,"uniform");
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
 save('output/feats4T.mat', 'features', '-mat' );

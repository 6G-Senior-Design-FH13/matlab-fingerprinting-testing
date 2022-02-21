TR = stlread("cornerWithFloor.stl");
viewer = siteviewer('SceneModel',TR);
mapFileName = "cornerWithFloor.stl";
tx = txsite("cartesian", ...
    "AntennaPosition",[0; 0; 2], ...
    "TransmitterFrequency",1e11);
show(tx,"ShowAntennaHeight",false) 

rx = rxsite("cartesian", ...
    "AntennaPosition",[4; 5; 1]);
show(rx,"ShowAntennaHeight",false)

pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","wood"); 

r = raytrace(tx,rx,pm);
r = r{1};
plot(r)

snrs = 10; 
chanBW = "CBW40"; 
txArraySize = [1 1];
rxArraySize = [1 1];

rays = raytrace(tx,rx,pm,"Map",mapFileName);

cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
    "NumTransmitAntennas",1, ...
    "SecureHELTF",false, ...
    "GuardInterval",1.6);
cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);

[features,labels] = dlPositioningGenerateDataSet(rays,rx,tx,cfg,snrs);

save('features_out.mat', 'features', '-mat' );
position=labels.position;
save('labels_out.mat', 'position', '-mat');
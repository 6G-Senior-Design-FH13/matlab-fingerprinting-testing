TR = stlread("office.stl");
mapFileName = "office.stl";
scale = 0.9;
scaledPts = TR.Points * scale;
TR_scaled = triangulation(TR.ConnectivityList,scaledPts);

viewer = siteviewer('SceneModel',TR_scaled);



tx = txsite("cartesian", ...
    "AntennaPosition",[0.01; 1; 2.5], ...
    "TransmitterFrequency",100e9); %100GHz max
show(tx,"ShowAntennaHeight",false) 

rx = rxsite("cartesian", ...
    "AntennaPosition",[4.5; 3.5; 0.8]);
show(rx,"ShowAntennaHeight",false)

txArraySize = [1 1];
rxArraySize = [1 1];

los(tx,rx);


pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","image", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","wood"); 
r = raytrace(tx,rx,pm);
% r = r{1};
plot(r{1});

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
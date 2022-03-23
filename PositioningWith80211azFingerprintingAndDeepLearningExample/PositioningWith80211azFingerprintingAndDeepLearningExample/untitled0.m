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

mapFileName = "office.stl";
TR = stlread("office.stl");
scale = 0.9;
scaledPts = TR.Points * scale;
TR_scaled = triangulation(TR.ConnectivityList,scaledPts);

viewer = siteviewer('SceneModel',TR_scaled);

xrx = [0.1 4.4];
yrx = [0.1 6.9];
zrx = 2.1;
antPosAP = [kron(xrx, ones(1, length(yrx))); ...
          repmat(yrx, 1, length(xrx)); ...
          zrx*ones(1, length(xrx)*length(yrx))];
snrs = 10; 
chanBW = "CBW40"; 
txArraySize = [1 1];
rxArraySize = [1 1];
fc = 25e9; % Set the carrier frequency (Hz)
lambda = physconst("lightspeed")/fc;

txArray = arrayConfig("Size", txArraySize, "ElementSpacing", 2*lambda);

tx = txsite("cartesian", ...
    "Antenna", txArray, ...
    "AntennaPosition", antPosAP,...
    "TransmitterFrequency", fc);
show(tx,"ShowAntennaHeight",false);
rx = rxsite("cartesian", ...
    "AntennaPosition",[1; 2; 2]);
show(rx,"ShowAntennaHeight",false)

los(tx(1),rx);
los(tx(2),rx);
los(tx(3),rx);
los(tx(4),rx);

pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","image", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","wood"); 
% r = raytrace(tx,rx,pm);
% plot(r{1});

rays1 = raytrace(tx(1),rx,pm,"Map",mapFileName);
rays2 = raytrace(tx(2),rx,pm,"Map",mapFileName);
rays3 = raytrace(tx(3),rx,pm,"Map",mapFileName);
rays4 = raytrace(tx(4),rx,pm,"Map",mapFileName);

cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
    "NumTransmitAntennas",1, ...
    "SecureHELTF",false, ...
    "GuardInterval",1.6);
cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);

% y = dlPositioningGenerateDataSet(rays,rx,tx,cfg,snrs);
[r1, y1] = dlPositioningGenerateDataSet(rays1,rx,tx(1),cfg,snrs);
[r2, y2] = dlPositioningGenerateDataSet(rays2,rx,tx(2),cfg,snrs);
[r3, y3] = dlPositioningGenerateDataSet(rays3,rx,tx(3),cfg,snrs);
[r4, y4] = dlPositioningGenerateDataSet(rays4,rx,tx(4),cfg,snrs);

A = [2*(tx(2).AntennaPosition(1,1) - tx(1).AntennaPosition(1,1)) 2*(tx(2).AntennaPosition(2,1) - tx(1).AntennaPosition(2,1)) 2*(tx(2).AntennaPosition(3,1) - tx(1).AntennaPosition(3,1));
    2*(tx(3).AntennaPosition(1,1) - tx(1).AntennaPosition(1,1)) 2*(tx(3).AntennaPosition(2,1) - tx(1).AntennaPosition(2,1)) 2*(tx(3).AntennaPosition(3,1) - tx(1).AntennaPosition(3,1));
    2*(tx(4).AntennaPosition(1,1) - tx(1).AntennaPosition(1,1)) 2*(tx(4).AntennaPosition(2,1) - tx(1).AntennaPosition(2,1)) 2*(tx(4).AntennaPosition(3,1) - tx(1).AntennaPosition(3,1))];
v = [tx(1).AntennaPosition(1,1)^2 + tx(1).AntennaPosition(2,1)^2 + tx(1).AntennaPosition(3,1)^2 - (tx(2).AntennaPosition(1,1)^2 + tx(2).AntennaPosition(2,1)^2 + tx(2).AntennaPosition(3,1)^2) - (r2^2 - r1^2);
    tx(1).AntennaPosition(1,1)^2 + tx(1).AntennaPosition(2,1)^2 + tx(1).AntennaPosition(3,1)^2 - (tx(3).AntennaPosition(1,1)^2 + tx(3).AntennaPosition(2,1)^2 + tx(3).AntennaPosition(3,1)^2) - (r3^2 - r1^2);
    tx(1).AntennaPosition(1,1)^2 + tx(1).AntennaPosition(2,1)^2 + tx(1).AntennaPosition(3,1)^2 - (tx(4).AntennaPosition(1,1)^2 + tx(4).AntennaPosition(2,1)^2 + tx(4).AntennaPosition(3,1)^2) - (r4^2 - r1^2)];
position = pinv(A'*A)*A'*v;
disp(position);
% disp(rx.AntennaPosition);
% figure(3);
% pspectrum(y1,40e6,'spectrogram', 'TimeResolution',0.000001, 'OverlapPercent',99);




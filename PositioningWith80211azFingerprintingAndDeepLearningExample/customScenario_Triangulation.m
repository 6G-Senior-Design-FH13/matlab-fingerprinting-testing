mapFileName = "office.stl";
TR = stlread("office.stl");
scale = 0.9;
scaledPts = TR.Points * scale;
TR_scaled = triangulation(TR.ConnectivityList,scaledPts);

viewer = siteviewer('SceneModel',TR_scaled);

xt = 2;
yt = 1; 
zt = 1;
xr = 1;
yr = 1;
zr = 1;

tx = txsite("cartesian", ...
    "AntennaPosition",[xt; yt; zt], ...
    "TransmitterFrequency",100e9); %100GHz max
show(tx,"ShowAntennaHeight",false) 

show(tx,"ShowAntennaHeight",false) 

rx = rxsite("cartesian", ...
    "AntennaPosition",[xr; yr; zr]);
show(rx,"ShowAntennaHeight",false)

los(tx,rx);

pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","image", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","wood"); 
r = raytrace(tx,rx,pm);
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

fname = ['f' '_' num2str(xt,'%d') '-' num2str(yt,'%d') '-' num2str(zt,'%d') '_' num2str(xr,'%d') '-' num2str(yr,'%d') '-' num2str(zr,'%d') '.mat'];
lname = ['l' '_' num2str(xt,'%d') '-' num2str(yt,'%d') '-' num2str(zt,'%d') '_' num2str(xr,'%d') '-' num2str(yr,'%d') '-' num2str(zr,'%d') '.mat'];


save(fname, 'features', '-mat' );
position=labels.position;
save(lname, 'position', '-mat');


% triangulation
a=rays{1,1}(1,1).PropagationDistance * cos(rays{1,1}(1,1).AngleOfDeparture(1)*pi/180) * cos(rays{1,1}(1,1).AngleOfDeparture(2)*pi/180);
b=rays{1,1}(1,1).PropagationDistance * sin(rays{1,1}(1,1).AngleOfDeparture(1)*pi/180) * cos(rays{1,1}(1,1).AngleOfDeparture(2)*pi/180);
c=rays{1,1}(1,1).PropagationDistance * sin(rays{1,1}(1,1).AngleOfDeparture(2)*pi/180);
array_ab = [a;b;c];
Rx_location = array_ab + rays{1,1}(1,1).TransmitterLocation;


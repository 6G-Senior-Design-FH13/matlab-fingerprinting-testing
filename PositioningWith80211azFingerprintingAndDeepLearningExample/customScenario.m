TR = stlread("corner.stl");
viewer = siteviewer('SceneModel',TR);

tx = txsite("cartesian", ...
    "AntennaPosition",[0; 0; 2], ...
    "TransmitterFrequency",2.8e9);
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
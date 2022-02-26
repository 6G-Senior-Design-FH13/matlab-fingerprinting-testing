TR = stlread("office.stl");
mapFileName = "office.stl";
scale = 0.9;
scaledPts = TR.Points * scale;
TR_scaled = triangulation(TR.ConnectivityList,scaledPts);

%viewer = siteviewer('SceneModel',TR_scaled);

xt = .1;
yt = .1;
zt = 2.5;
%xr = 2; %.1 to 4.5 45 samps
%yr = 7.1; %.1 to 7.1 71 samps
%zr = 2.5; %1 to 2.5 16 samps total 51120
count = 0;
feats = zeros(48,51120);
labels = zeros(3,51120);
for xr=.1:.1:4.5
    for yr=.1:.1:7.1
        for zr=1:.1:2.5
            count=count+1;
            
            tx = txsite("cartesian", ...
                "AntennaPosition",[xt; yt; zt], ...
                "TransmitterFrequency",100e9); %100GHz max
            %show(tx,"ShowAntennaHeight",false) 
            
            rx = rxsite("cartesian", ...
                "AntennaPosition",[xr; yr; zr]);
            %show(rx,"ShowAntennaHeight",false)
            
            %los(tx,rx);
            
            pm = propagationModel("raytracing", ...
                "CoordinateSystem","cartesian", ...
                "Method","image", ...
                "MaxNumReflections",2, ...
                "SurfaceMaterial","wood"); 
            
            r = raytrace(tx,rx,pm);
            if(~isempty(r))
                % r = r{1};
                %plot(r{1});
                
                snrs = 10; 
                chanBW = "CBW40"; 
                txArraySize = [1 1];
                rxArraySize = [1 1];
                
                rays = raytrace(tx,rx,pm,"Map",mapFileName);
                rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
                cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
                    "NumTransmitAntennas",1, ...
                    "SecureHELTF",false, ...
                    "GuardInterval",1.6);
                cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);
                
                [features] = dlPositioningGenerateDataSet(rays,rx,tx,cfg,snrs);
                
                feats(:,count) = features;
                labels(:,count)= [xr yr zr];
                %fname = ['output/f' '_' num2str(xr,'%.1f') '-' num2str(yr,'%.1f') '-' num2str(zr,'%.1f') '.mat'];
                %lname = ['output/l' '_' num2str(xr,'%.1f') '-' num2str(yr,'%.1f') '-' num2str(zr,'%.1f') '.mat'];
                
                
                %save(fname, 'features', '-mat' );
                %position=labels.position;
                %save(lname, 'position', '-mat');
            end   
        end
    end
    save('output/feats.mat', 'feats', '-mat' );
    save('output/labels.mat', 'labels', '-mat' );
end


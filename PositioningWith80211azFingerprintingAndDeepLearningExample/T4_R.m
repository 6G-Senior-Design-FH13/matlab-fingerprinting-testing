    mapFileName = "office.stl";
    viewer = siteviewer("SceneModel",mapFileName,"Transparency",0.25);
    distribution = "uniform";
    txArraySize = [4 1]; % Linear transmit array
    rxArraySize = [1 1]; % Linear receive array
    chanBW = "CBW40"; 
    staSeparation = .5; % STA separation, in meters, used only when the distribution is uniform
    numSTAs = 480;  
    S = RandStream("mt19937ar","Seed",5489); % Set the RNG for reproducibility.
    RandStream.setGlobalStream(S);
    if distribution == "uniform"
        [APs,STAs] = DEFdlPositioningCreateEnvironment(txArraySize,rxArraySize,staSeparation,"uniform");
    else
        [APs,STAs] = DEFdlPositioningCreateEnvironment(txArraySize,rxArraySize,numSTAs,"random");
    end
    pm = propagationModel("raytracing", ...
        "CoordinateSystem","cartesian", ...
        "SurfaceMaterial","wood", ...
        "MaxNumReflections",2);
    %rays = raytrace(APs,STAs,pm,"Map",mapFileName);
    snr = 10; 
    show(APs)
    show(STAs(30),'IconSize',[32 32]);
    plot([rays{:,30}],'ColorLimits',[50 95]);
    cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
        "NumTransmitAntennas",prod(txArraySize), ...
        "SecureHELTF",false, ...
        "GuardInterval",1.6);
    cfg.User{1}.NumSpaceTimeStreams = prod(txArraySize);
    
    [features,labels] = DEFdlPositioningGenerateDataSet(rays,STAs,APs,cfg,snr);
    
    lp = labels.position;
    
    if distribution == "uniform"
        save('output/feats4T_.5R_3refl.mat', 'features', '-mat' );
        save('output/labels4T_.5R_3refl.mat', 'lp', '-mat' );
    else
        save('output/feats4T_480R_2refl.mat', 'features', '-mat' );
        save('output/labels4T_480R_2refl.mat', 'lp', '-mat' );
    end


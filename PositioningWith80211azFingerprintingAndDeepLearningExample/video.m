mapFileName = "office.stl";
    viewer = siteviewer("SceneModel",mapFileName,"Transparency",0.25);
    daspect([1,1,.3]);axis tight;
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    distribution = "uniform";
    txArraySize = [4 1]; % Linear transmit array
    rxArraySize = [1 1]; % Linear receive array
    chanBW = "CBW40"; 
    staSeparation = 1; % STA separation, in meters, used only when the distribution is uniform
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
    rays = raytrace(APs,STAs,pm,"Map",mapFileName);
    snr = 10; 
    show(APs)
    show(STAs(30),'IconSize',[32 32]);
    plot([rays{:,30}],'ColorLimits',[50 95]);
    

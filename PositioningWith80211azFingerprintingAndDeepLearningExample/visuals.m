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
    cfg = heRangingConfig('ChannelBandwidth',chanBW, ...
        "NumTransmitAntennas",prod(txArraySize), ...
        "SecureHELTF",false, ...
        "GuardInterval",1.6);
    txWaveform=cfg.ChannelBandwidth;
    snr = 10; 
    show(APs)
    show(STAs,'ShowAntennaHeight',false,'IconSize',[16 16]);
    figure()
    array=features(:,1,2,1);
    plot(array);
    hold on
    r = randi(480,10,1);
    for i=1:10
        array=features(:,1,2,r(i));
        plot(array);
    end
    hold off
    title('CIRs of 10 Random Receivers')
    set(gca, 'YScale', 'log')
    xlabel('Time Slice "n"') 
    ylabel('log(Channel Impulse Response)') 
TR = stlread("office.stl");

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


% 5G waveform generate

carrier = nrSCSCarrierConfig('NSizeGrid',100);
bwp = nrWavegenBWPConfig('NStartBWP',carrier.NStartGrid+10);
ssb = nrWavegenSSBurstConfig('BlockPattern','Case A');
pdcch = nrWavegenPDCCHConfig('AggregationLevel',2,'AllocatedCandidate',4);
coreset = nrCORESETConfig;
coreset.FrequencyResources = [1 1 1 1];
coreset.Duration = 3;
ss = nrSearchSpaceConfig;
ss.NumCandidates = [8 4 0 0 0];
pdsch = nrWavegenPDSCHConfig( ...
    'Modulation','16QAM','TargetCodeRate',658/1024,'EnablePTRS',true);
dmrs = nrPDSCHDMRSConfig('DMRSTypeAPosition',3);
pdsch.DMRS = dmrs;
ptrs = nrPDSCHPTRSConfig('TimeDensity',2);
pdsch.PTRS = ptrs;
csirs = nrWavegenCSIRSConfig('RowNumber',4,'RBOffset',10);
cfgDL = nrDLCarrierConfig( ...
    'FrequencyRange','FR1', ...
    'ChannelBandwidth',40, ...
    'NumSubframes',20, ...
    'SCSCarriers',{carrier}, ...
    'BandwidthParts',{bwp}, ...
    'SSBurst',ssb, ...
    'CORESET',{coreset}, ...
    'SearchSpaces',{ss}, ...
    'PDCCH',{pdcch}, ...
    'PDSCH',{pdsch}, ...
    'CSIRS',{csirs});

txWaveform = nrWaveformGenerator(cfgDL);

% Reconstruct Channel Impulse Response Using CDL Channel Path Filters (google)

v = 15.0;                    % UE velocity in km/h
fc = 4e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
 
cdl = nrCDLChannel;
cdl.DelayProfile = 'CDL-D';
cdl.DelaySpread = 10e-9;
cdl.CarrierFrequency = fc;
cdl.MaximumDopplerShift = fd;

cdl.TransmitAntennaArray.Size = [txArraySize 1 1 1];
cdl.ReceiveAntennaArray.Size = [rxArraySize 1 1 1];

SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumTransmitAntennas;
 
% txWaveform = complex(randn(T,Nt),randn(T,Nt));

[rxWaveform,pathGains] = cdl(txWaveform);

pathFilters = getPathFilters(cdl);

[offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
%  mag = magnitude of the channel impulse response
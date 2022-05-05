%JC 3/7/08
%UWB Transmitter/Receiver using BPSK Continuous Wave
%Run from editor debug(F5).
%m-file uses random data which BPSK modulates a
%carrier to construct a  BPSK UWB transmitter. The receiver demodulates the BPSK UWB
%signal (assumming perfect sync) and the filtered data is recovered.
%NOTE:You should assume that the variables are multiplied by ONE MILLION to theoritically
%give values such as fcarr=4.0GHz,N=1.3GHz and so on. The numbers shown here, when
%upscaled, signify pulse widths of ~ 750 pico seconds (1/1.3GHz). In other words, what
%I'm trying to say is you have a system that can operate in the 3.1-10.6 GHz
%UWB band(at data rates up to 1.3 Gbps) by changing fcarr and the bandwidth can be changed
%by varing the time of the data rate or N. It would be nice to do this in real time with 
%multipath, FEC, PLL's etc, but we do with what we have. Another thought would be to use 
%a QPSK format to increase the data rate. I have provided references at the end of this
%file. They should give you a good idea of this particuliar technique and some possible
%differences between it and impulse response UWB and MB-OFDM UWB. The program is somewhat
%rudimentary but I hope it is clear enough to be helpful.

clear;
fcarr=4e3;              % Carrier frequency
N = 1300;		        % Number of symbols (data rate)
fs = 15.6*1e3;		    % Sampling frequency
Fn = fs/2;              % Nyquist frequency
Ts = 1/fs;	            % Sampling time = 1/fs
T = 1/N;		        % Symbol time
randn('state',0);       % Keeps transmitted data from changing on reruns
t = [0:Ts:(1-Ts)];      % Time vector
%==========Transmitter================================================
symbols = sign(randn(N,1))';%generate N random binary symbols(data), +/- 1
%(notice transpose')
symbols1 = ones(T/Ts,1)*symbols;
s2 = symbols1(:); 	%data @ +/-1

%generate carrier wave
%sine wave
%2 pi fc t is written as below
twopi_fc_t=(1:fs)*2*pi*fcarr/fs; 
a=2;  %set amplitude
f_c = a * sin(twopi_fc_t); %unmodulated carrier frequency
%======================================================================
s2=s2';            %transpose for matrix match

f_c1=f_c.*s2;     %multiply carrier and +/- signal data to get modulated BPSK carrier
%At this point in the transmitter you would add a bandpass filter to meet
%FCC UWB spectral mask of -41.25 dBm/MHz EIRP.

%========================================================================
%take FFT of modulated carrier(f_c1) [spectrum analyzer]
y=f_c1;
NFFY=2.^(ceil(log(length(y))/log(2)));
FFTY=fft(y,NFFY);%pad with zeros
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(y);
f1=(0:NumUniquePts-1)*2*Fn/NFFY;

%plot frequency domain
subplot(3,2,4); plot(f1,MY);xlabel('');ylabel('AMPLITUDE');
axis([0 8000 -.5 .5]);%zoom in/out
title('FREQUENCY DOMAIN PLOTS');
grid on;
subplot(3,2,6); plot(f1,20*log10(abs(MY).^2));xlabel('FREQUENCY(MHz)');ylabel('DB');
axis([0 8000 -100 -20]);
grid on;
title(' UNFILTERED BPSK MODULATED UWB SIGNAL ')

%plot data [ocilloscope-time domain]
figure(1)
subplot(3,2,1)
plot(t,s2)
axis([0 4.5e-3 -1.2 1.2])
grid on
title('TRANSMITTED DATA')

subplot(3,2,3) 
plot(t,f_c1)
axis([0 4.5e-3 -4 4])
grid on
title('MODULATOR OUT 0/180 DEGREE PHASE SHIFTS ')

subplot(3,2,5) 
plot(t,f_c)
axis([0 4.5e-3 -4 4])
grid on
title('CARRIER FREQUENCY')

%===============Receiver===================================================
%Assume perfect sync
data1=f_c.*f_c1; %multiply carrier with modulated output
%At this point in the receiver, you would add a low pass filter and use a
%one to three bit flash ADC as the signal level will be quiet small 
%(depending on distance)even with a low noise pre amp.
subplot(3,2,2)
plot(t,data1)
axis([0 4.5e-3 -5.2 5.2])
grid on
title('UNFILTERED RECEIVED DATA')

%References
%http://www.mwee.com/590028720 "All Ultra-Wideband (UWB) systems are not
%created equal"
%http://electronicdesign.com   "Third UWB Method Solves The Home Network
%Problem"
%http://www.wirelessnetdesignline.com  "Comprehensive UWB product testing"


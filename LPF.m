clc
clear all
close all

fs=350; %sampling frequency
nyquist=fs/2; %nyquist frequency
endpoint=7; %number of time points on time axis
timevec=(0:fs*endpoint-1)/fs; %time axis 
points= length(timevec); %number of time points in the signal
ybrown=cumsum(randn(points,1)); %generate signal with brownian noise
ywhite=randn(points,1); %generate signal with white noise
yspike=sin(2*pi*50*timevec); %generate a large spike at 50 Hz
y=ybrown+ywhite+yspike; %signal on which low pass filter is to be applied
                        %(Original Signal)
subplot(231)
plot(timevec,ybrown,'k') 
xlabel('Time (s)')
title('Signal with Brownian Noise') 
subplot(232)
plot(timevec,ywhite,'k'); 
xlabel('Time (s)')
title('Signal with White Noise') 

freqvec=linspace(0,nyquist,floor(points/2)+1); %frequency axis 
yspikef=abs(fft(yspike)/fs); %large spike in frequency domain
subplot(233)
plot(freqvec,yspikef(1:length(freqvec)),'k') 
xlabel('Frequency (Hz)')
title('Spike at 50Hz')

subplot(234)
plot(timevec,y,'b')
xlabel('Time (s)')
title('Original Signal') 

yspectrum=abs(fft(y)/points).^2; %power spectrum of original signal
subplot(235)
plot(freqvec,yspectrum(1:length(freqvec)),'b')
set(gca,'yscale','log')
xlabel('Frequency (Hz)')
title('Power Spectrum of Original Signal')

%Apply Low Pass Filter (firls)
cutoff= 30; %cutoff frequency (Hz)
transw=20/100; %transition width (percentage)
order=round(7*fs/cutoff); %order of the filter * sampling frequency

shape=[1 1 0 0]; %amplitude vector for firls
freq= [0 cutoff cutoff+cutoff*transw nyquist]/(nyquist); %frequency vector for firls
fkern=firls(order,freq,shape); %genrate filter kernel by firls
fkernspectrum= abs(fft(fkern,points)).^2; %spectrum of filter kernel

figure
subplot(221)
plot((-order/2:order/2)/fs,fkern,'k') %length of the kernel = order of the filter
xlabel('Time (s)')
title('Filter Kernel')

subplot(222)
plot(freq*nyquist,shape,'r') %ideal low pass filter kernel
hold on
plot(freqvec,fkernspectrum(1:length(freqvec)),'k') %real low pass filter kernel
xlabel('Frequency (Hz)')
title('Filter kernel spectrum')
legend('Ideal','Real')

yfilt=filtfilt(fkern,1,y); %apply filter to original signal

subplot(223)
plot(timevec,y,'b') %original signal
hold on
plot(timevec,yfilt,'m') %filtered signal
xlabel('Time (s)')
title('Filtered Signal')
legend('Original','Filtered')

yfiltspectrum=abs(fft(yfilt)/points).^2; %power spectrum of filtered signal

subplot(224)
plot(freqvec,yspectrum(1:length(freqvec)),'b') %original signal
hold on
plot(freqvec,yfiltspectrum(1:length(freqvec)),'m') %filtered signal
set(gca,'yscale','log')
xlabel('Frequency (Hz)')
title('Power Spectrum of Filtered Signal')
legend('Original','Filtered')



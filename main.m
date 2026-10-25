% Fourier and amplitude modulation experiment, reproducible without audio files.
clear; close all;
if ~exist('output','dir'), mkdir('output'); end
fs=96000; N=16384; t=(0:N-1)/fs;
f=(-N/2:N/2-1)*fs/N;
figure('visible','off');
widths=[0.001,0.003,0.006];
for k=1:length(widths)
 pulse=double(abs(t-0.05)<widths(k)/2);
 spectrum=fftshift(fft(pulse))/fs;
 subplot(3,2,2*k-1);plot(t,pulse);xlim([.04,.06]);ylabel('Amplitude');title(sprintf('Width %.3f s',widths(k)));
 subplot(3,2,2*k);plot(f,abs(spectrum));xlim([-5000,5000]);ylabel('|X(f)|');xlabel('Hz');
end
print('output/pulse_spectra.png','-dpng');
x=[1,2,3,4];h=[.25,.5,.25];L=length(x)+length(h)-1;
direct=conv(x,h);frequency=real(ifft(fft(x,L).*fft(h,L)));
assert(max(abs(direct-frequency))<1e-12);
message=.55*sin(2*pi*440*t)+.25*sin(2*pi*1000*t)+.1*sin(2*pi*2400*t);
carrier=12000;modulated=(1+message).*cos(2*pi*carrier*t);
% Coherent demodulation and zero-phase ideal FFT low-pass reference.
mixed=2*modulated.*cos(2*pi*carrier*t);
mask=abs(f)<3500;
recovered=real(ifft(ifftshift(fftshift(fft(mixed)).*mask)))-1;
% Edge leakage is excluded from the central-window metric.
region=1000:N-1000;rmse=sqrt(mean((message(region)-recovered(region)).^2));
assert(rmse<.05);
figure('visible','off');
subplot(3,1,1);plot(t,message);xlim([0,.015]);title('Message');ylabel('Amplitude');
subplot(3,1,2);plot(f,abs(fftshift(fft(modulated))));xlim([-20000,20000]);title('AM spectrum');xlabel('Hz');
subplot(3,1,3);plot(t,message,t,recovered);xlim([0,.015]);legend('Original','Recovered');xlabel('Seconds');
print('output/am_reconstruction.png','-dpng');
fid=fopen('output/metrics.txt','w');fprintf(fid,'Convolution max error: %.12g\nReconstruction RMSE: %.8f\n',max(abs(direct-frequency)),rmse);fclose(fid);
fprintf('PASS convolution and AM reconstruction; RMSE %.8f\n',rmse);

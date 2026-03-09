clear

% RF bandwidth (Hz)
bandwidth = 2000;

% number of side lobes (sharpness) 
nlobes = 1;

% main lobe duration (us)
duration = 2e6/bandwidth;

% number of points (us)
np = (1+nlobes) * duration;
np = 2*ceil(np/2); % even integer

% time vector (us)
t = -np/2:np/2; % center *exactly* 0
dt = t(2)-t(1); % us

% RF shape (windowed sinc)
a = sinc(bandwidth*t*1e-6);
a = a.*(0.54+0.46*cospi(2*t/np));

% hard pulse RF
%a=double(cospi(2*t/np)>0);

% calculate main lobe width
hi = np/2;
lo = np/2;
for k = 1:np/2
    if a(hi)>0; hi = hi+1; end
    if a(lo)>0; lo = lo-1; end
end
if hi<np; hi = interp1(a([hi-1 hi]),[hi-1 hi],0); end
if lo> 0; lo = interp1(a([lo lo+1]),[lo lo+1],0); end
WIDTH = (hi-lo) * dt; 

% half pulse RF
%for k = 1+np/2:np
%     a(k)=0; 
%end

%% Fourier domain 
N = 10*np; % pad to N points and extract center
A = fft(circshift([a zeros(1,N-np)],-np/2));
A = circshift(A,np-1);
A = A(1:2*np);
A = A / max(A);

% sampling frequency
Fs = (1e6/dt)/N; % Hz
f = (1-numel(A)/2:numel(A)/2) * Fs;

% calculate fwhm ~ 2/WIDTH
hi = numel(A)/2;
lo = numel(A)/2;
for k = 1:numel(A)-1
    if real(A(hi+1))>0.5; hi = hi+1; end
    if real(A(lo-1))>0.5; lo = lo-1; end       
    if real(A(hi+1))<=0.5 && real(A(lo-1))<=0.5; break; end
end
hi = interp1(real(A(hi+(0:1))),hi+(0:1),0.5);
lo = interp1(real(A(lo-(0:1))),lo-(0:1),0.5);
FWHM = (hi-lo) * Fs;

%% display (RF pulse, freq response, fid decay)
subplot(1,2,1)
cplot(t,a);
xlabel('Time (us)');
title('RF waveform')
grid on
xlim([min(t) max(t)]); yticks(0.25*(-4:4));
for k=[-1 1]
    line([k k]*WIDTH/2,ylim,'linestyle',':','color','black');
end
text(-450,-0.05,['MAIN LOBE WIDTH = ' num2str(WIDTH*1e-3,'%.2f') 'ms'])
legend({'real','imag'});

subplot(1,2,2)
cplot(f,A);
xlabel('Frequency (Hz)');
title('Response');
grid on
xlim([-4e6 4e6]/WIDTH); yticks(0.25*(-4:4));
for k=[-1 1]
    line([k k]*FWHM/2,ylim,'linestyle',':','color','black');
end
text(1100,0.525,['FWHM = ' num2str(FWHM*1e-3,'%.2f') 'kHz'])


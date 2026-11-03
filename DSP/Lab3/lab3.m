%%
% z-transform 
num=[1,2,1.3,0,4,1];
den=[0.65,0,2.3,3,2,1];
[z,p,k]=tf2zp(num,den);
zplane(z,p)

%%
%z-transform 
z =[
  -1.4857 + 1.0130i
  -1.4857 - 1.0130i
   0.6090 + 0.9394i
   0.6090 - 0.9394i
  -0.2467 + 0.0000i
]; % zeros 
p =[
   0.6212 + 1.9665i
   0.6212 - 1.9665i
  -0.7748 + 0.0000i
  -0.2338 + 0.6420i
  -0.2338 - 0.6420i
]; % poles 

k =1.5385; % gain 

[num,den]=zp2tf(z,p,k);

% Normalize
num = num / den(1);
den = den / den(1);

% Display H(z)
Hz = tf(num, den, -1)


%%
% frequency response 
num=[0.008,-0.033,0.05,-0.033,0.008]
den=[1,2.37,2.7,1.6,0.4]
w=20 % range of frequency 
x=freqz(num,den,w)
subplot(211)
plot(abs(x))
grid on
subplot(212)
plot(angle(x))
grid on 

%%
%FIR filter design 
N=input('enter the length of the filter')
wc=input('cutoff frequency=')
h=fir1(N,wc/pi,kaiser(N+1)) % response of filter
freqz(h)

%change blackman with rectwin, hanning etc

%%
% FIR filter design from kaiser window 
N = input('Enter filter length: ');
wc = input('Cutoff frequency: ');

beta = 5; % typical value (~60 dB attenuation)

h = fir1(N, wc/pi, kaiser(N+1, beta));

freqz(h)

%%
% IIR filter design 
rp=0.4 % passband ripple
rs=30 %stopband attenuation 
wp=0.2*pi %passband frequency 
ws=0.4*pi %stopband frequency 
[N,wn]=buttord(wp/pi , ws/pi, rp, rs) %finding order of the filter  
% [b,a]=butter(N,wn) %low pass 
[c,d]=butter(N,wn, 'high') %high pass
% freqz(b,a)
freqz(c,d)


%%
%FFT 
% Define the input sequence
x = [1.2,3.6,2.4,-2,1.6,2.1,0,0];

% Compute the DFT using the built-in fft function
X = fft(x,8);
subplot(411)
stem(real(X))
subplot(412)
stem(imag(X))
subplot(413)
stem(abs(X))
subplot(414)
stem(angle(X))

% Display the result
disp('The DFT of x[n] is:');
disp(X);

% To see the Magnitude and Phase 
magnitude = abs(X);
phase = angle(X);

fprintf('Magnitudes: [%.2f, %.2f, %.2f, %.2f]\n', magnitude);
fprintf('Phases (rad): [%.2f, %.2f, %.2f, %.2f]\n', phase);

%%
% FFT 
% Define the input sequence
x = [0, 2, 1, 3];

% Compute the DFT using the built-in fft function
X = fft(x);
subplot(411)
stem(real(X))
subplot(412)
stem(imag(X))
subplot(413)
stem(abs(X))
subplot(414)
stem(angle(X))

% Display the result
disp('The DFT of x[n] is:');
disp(X);

% To see the Magnitude and Phase 
magnitude = abs(X);
phase = angle(X);

fprintf('Magnitudes: [%.2f, %.2f, %.2f, %.2f]\n', magnitude);
fprintf('Phases (rad): [%.2f, %.2f, %.2f, %.2f]\n', phase);


%%
%%
% FFT ANALYSIS
clc; clear; close all;

%% 🔹 Case 1
x1 = [1.2, 3.6, 2.4, -2, 1.6, 2.1, 0, 0];
N1 = length(x1);

X1 = fft(x1, N1);

subplot(4,1,1)
stem(0:N1-1, real(X1), 'filled', 'LineWidth', 1.5)
title('Real Part'); xlabel('k'); ylabel('Re{X(k)}'); grid on;

subplot(4,1,2)
stem(0:N1-1, imag(X1), 'filled', 'LineWidth', 1.5)
title('Imaginary Part'); xlabel('k'); ylabel('Im{X(k)}'); grid on;

subplot(4,1,3)
stem(0:N1-1, abs(X1), 'filled', 'LineWidth', 1.5)
title('Magnitude'); xlabel('k'); ylabel('|X(k)|'); grid on;

subplot(4,1,4)
stem(0:N1-1, angle(X1), 'filled', 'LineWidth', 1.5)
title('Phase'); xlabel('k'); ylabel('Angle (rad)'); grid on;

disp('Input x[n] = '); disp(x1)

fprintf('\nDFT X[k]:\n');
disp(X1.')

fprintf('\n k    Magnitude    Phase(rad)\n');
for k = 1:N1
    fprintf('%2d    %8.3f      %8.3f\n', k-1, abs(X1(k)), angle(X1(k)));
end


%% 🔹 Case 2
x2 = [0, 2, 1, 3];
N2 = length(x2);

X2 = fft(x2, N2);

subplot(4,1,1)
stem(0:N2-1, real(X2), 'filled', 'LineWidth', 1.5)
title('Real Part'); xlabel('k'); ylabel('Re{X(k)}'); grid on;

subplot(4,1,2)
stem(0:N2-1, imag(X2), 'filled', 'LineWidth', 1.5)
title('Imaginary Part'); xlabel('k'); ylabel('Im{X(k)}'); grid on;

subplot(4,1,3)
stem(0:N2-1, abs(X2), 'filled', 'LineWidth', 1.5)
title('Magnitude'); xlabel('k'); ylabel('|X(k)|'); grid on;

subplot(4,1,4)
stem(0:N2-1, angle(X2), 'filled', 'LineWidth', 1.5)
title('Phase'); xlabel('k'); ylabel('Angle (rad)'); grid on;

disp('Input x[n] = '); disp(x2)

fprintf('\nDFT X[k]:\n');
disp(X2.')

fprintf('\n k    Magnitude    Phase(rad)\n');
for k = 1:N2
    fprintf('%2d    %8.3f      %8.3f\n', k-1, abs(X2(k)), angle(X2(k)));
end

% Discrete-Time Unit Impulse Signal 
% hold on 
% for(n=-5:5)
%     if(n==0)
%         stem(n,1,'r')
%     else 
%         stem(n,0,'r')
%     end
% end
% hold off 
% xlabel('n')
% ylabel('x[n]')
% title('Discrete-time unit impulse signal')

% Continous-Time Unit Impulse Signal
% t = -5:0.01:5; 
% plot(t, zeros(size(t)), 'r', "LineWidth", 2.5)
% hold on 
% line([0, 0], [0, 1], 'color', 'r', "LineWidth", 2.5)
% hold off 
% xlabel('t')
% ylabel('x(t)')
% title('Continous-time unit impulse signal')

% Continous-Time Unit Step Signal
% t = -5:0.01:5; 
% x = double(t>=0); 
% plot(t, x, 'r', "LineWidth", 2.5)
% xlabel('t')
% ylabel('x(t)')
% title('Continous-time unit step signal')

% Discrete-Time Unit Step Signal
% hold on 
% for(n=-5:5)
%     if(n<0)
%         stem(n, 0)
%     else
%         stem(n, 1, "LineWidth", 2.5)
%     end
% end
% hold off 
% xlabel('n')
% ylabel('x[n]')
% title('Discrete-time unit step signal')

% Continous-Time Ramp signal 
% m=1 and c=2 
% t = -5:0.01:5; 
% m = input('m=');
% c = input('c=');
% y = m*t + c; 
% plot(t, y, "LineWidth", 2.5)
% xlabel('t')
% ylabel('x(t)')
% title('Continous-time ramp signal')

% Discrete-Time Ramp Signal
% m = 2 and c = -1
% n = -5:5; 
% m = input('m=');
% c = input('c=');
% y = m*n + c; 
% stem(n, y, "LineWidth", 2.5)
% xlabel('n')
% ylabel('x[n]')
% title('Discrete-time ramp signal')

% Continous-Time Exponential Signal 
% a = 1 and c = 2 
% a = -1 and c = 2
% t = -5:0.01:5; 
% a = input('a=');
% c = input('c=');
% y = c*exp(a*t); 
% plot(t, y, "LineWidth", 2.5)
% xlabel('t')
% ylabel('x(t)')
% title('Continous-time exponential signal')

%Discrete-Time exponential Signal
% a = 0.3 and c = 1
% a = -0.5 and c = 1
% n = -5:5; 
% c = input('c=');
% a = input('a=');
% y = c*exp(a*n); 
% stem(n, y, "LineWidth", 2.5)
% xlabel('n')
% ylabel('x[n]')
% title('Discrete-time exponential signal')

% Continous-Time Sinusoidal signal
% hold on 
% t = -5:0.01:5;
% x = sin(t);
% y = cos(t);
% plot(t, x, "LineWidth", 2.5)
% plot(t, y, "LineWidth", 2.5)
% xlabel('t')
% ylabel('x(t) and y(t)')
% title('Continous-time sinusoidal signal')
% hold off

% Discrete-Time Sinusoidal signal 
% n = -5:1:5;
% x = sin(n);
% stem(n, x, "LineWidth", 2.5)
% xlabel('n')
% ylabel('x[n]')
% title('Discrete-time exponential signal')

% x = [1, 0, 1, 1];
% h = [1, -1, 0, 1];
% y = conv(x, h);
% subplot(3, 1, 1)
% stem(x,'filled', "LineWidth", 2.5)
% title('Input signal x[n]')
% xlabel('n')
% ylabel('x[n]')
% grid on 
% subplot(3, 1, 2)
% stem(h, 'filled', "LineWidth", 2.5)
% title('Impulse response h[n]')
% xlabel('n')
% ylabel('h[n]')
% grid on 
% subplot(3, 1, 3)
% stem(y, 'filled', "LineWidth", 2.5)
% title('Output signal y[n]')
% xlabel('n')
% ylabel('y[n]')
% grid on 


n = 0:10;
x = 0.5.^n;
h = ones(1, 11);
y = conv(x, h);
subplot(3, 1, 1)
stem(x,'filled', "LineWidth", 2.5)
title('Input signal x[n]')
xlabel('n')
ylabel('x[n]')
grid on 
subplot(3, 1, 2)
stem(h, 'filled', "LineWidth", 2.5)
title('Impulse response h[n]')
xlabel('n')
ylabel('h[n]')
grid on 
subplot(3, 1, 3)
stem(y, 'filled', "LineWidth", 2.5)
title('Output signal y[n]')
xlabel('n')
ylabel('y[n]')
grid on 
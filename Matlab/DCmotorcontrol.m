% N_G = [0, 0, 0, 1];
% D_G = [0.5, 100.5, 100, 0];
% G = tf(N_G, D_G);
% N_H = [1, 0];
% D_H = [0, 1]; 
% H = tf(N_H, D_H);
% N_L = [0, 0, 1, 0];
% D_L = [0.5, 100.5, 100, 0];
% L = tf(N_L, D_L);
% N_T = [0, 0, 0, 1];
% D_T = [0.5, 100.5, 101, 0];
% T = tf(N_T, D_T);
% 
% figure 
% bodeplot(T); 
% grid 
% 
% figure 
% nxquistplot(T);
% grid 
% 
% figure 
% w = logspace(-1, 2, 1000);
% [mag, phase] = bode(T, w);
% mag = squeeze(mag);
% phase = squeeze(phase);
% polarplot(deg2rad(phase), mag, 'LineWidth', 3);
% grid on 
% 
% figure 
% rlocusplot(T)
% grid 
% 
% 
% 
% 
% 




R1 = 10;
R2 = 20; 
C1 = 10e-6;
C2 = 200e-6; 
A = R1*R2*C1*C2; 
B = R1*C1 + R2*C2 + R1*C2; 
C = 1;
D = [A, B, C];
N = [0, 0, 1];
T = tf(N, D);

figure 
bodemag(T)
grid on 

figure
step(T)
grid on 

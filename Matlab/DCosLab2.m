% clear;
% DT = 2e-6; 
% order = 1; 
% ns = [0, 1];
% ds = [.03, 10];
% as = tf(ns, ds); 
% zbi = c2d(as, DT, 'tustin');
% zfh = c2d(as, DT, 'foh');
% zzh = c2d(as, DT, 'zoh');
% 
% figure
% step(as, 'r-', zzh, 'k-*', zfh, 'b--', zbi, 'c-.')
% grid;
% 
% figure
% bodemag(as, 'r-', zzh, 'k-*', zfh, 'b--', zbi, 'c-.')
% grid;


clear;
DT = 1; 
order = 1; 
ns = [0, 0, 1];
ds = [10.9, 1, 1.42];
as = tf(ns, ds); 
zbi = c2d(as, DT, 'tustin');
zfh = c2d(as, DT, 'foh');
zzh = c2d(as, DT, 'zoh');

figure
step(as, 'r-', zzh, 'k-*', zfh, 'b--', zbi, 'c-.')
grid;

figure
bodemag(as, 'r-', zzh, 'k-*', zfh, 'b--', zbi, 'c-.')
grid;

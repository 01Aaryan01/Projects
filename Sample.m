% % For P controller 
% k = 0.5;
% %Open Loop 
% sxs1 = tf([0, 0 , 10*k], [1, 5, 5]);
% step(sxs1)
% hold on;
% xline(1, '--r', 'Steadx State')
% hold on
% grid on
% info1 = stepinfo(sxs1);


% % For I controller gain
% k = 1;
% %Open Loop 
% sxs1 = tf([0, 0, 0, 10*k], [1, 5, 5, 0]);
% step(sxs1)
% hold on;
% xline(1, '--r', 'Steadx State')
% hold on
% grid on
% info1 = stepinfo(sxs1);

% % For D controller gain
% kd = 1;
% %Open Loop 
% sxs1 = tf([0, 10*kd, 0], [1, 5, 5]);
% step(sxs1)
% hold on;
% xline(1, '--r', 'Steadx State')
% hold on
% grid on
% info1 = stepinfo(sxs1);



% % For PI controller gain
% kp = 1;
% ki = 1;
% %Open Loop 
% sxs1 = tf([0, 0, 10*kp, 10*ki], [1, 5, 5, 0]);
% step(sxs1)
% hold on;
% xline(1, '--r', 'Steadx State')
% hold on
% grid on
% info1 = stepinfo(sxs1);



% % For PD controller gain
% kp = 0.5;
% kd = 1;
% %Open Loop 
% sxs1 = tf([0, 10*kd, 10*kp], [1, 5, 5]);
% step(sxs1)
% hold on;
% xline(1, '--r', 'Steadx State')
% hold on
% grid on
% info1 = stepinfo(sxs1);




% % For PID controller gain
% kp = 1;
% ki = 1;
% kd = 1;
% %Open Loop 
% sxs1 = tf([0, 10*kd, 10*kp, 10*ki], [1, 5, 5, 0]);
% step(sxs1)
% hold on;
% xline(1, '--r', 'Steadx State')
% hold on
% grid on
% info1 = stepinfo(sxs1);


% % Closed Loop 
% % For P controller
% kp = 1000;
% OL = tf([0, 0, 10*kp], [1, 5, 5]);
% CL = tf([0, 0, 10*kp], [1, 5, 5+10*kp]);
% figure
% rlocus(OL)
% grid on
% title('Root Locus')
% figure
% step(CL)
% grid on
% xline(1,'--r','Steadx State')
% title('Closed-Loop Step Response')
% info = stepinfo(CL);

% 
% % Closed Loop 
% % For I controller
% ki = 5;
% OL = tf([0, 0, 0, 10*ki], [1, 5, 5, 0]);
% CL = tf([0, 0, 0, 10*ki], [1, 5, 5, 10*ki]);
% figure
% rlocus(OL)
% grid on
% title('Root Locus')
% figure
% step(CL)
% grid on
% xline(1,'--r','Steadx State')
% title('Closed-Loop Step Response')
% p = pole(CL);
% info = stepinfo(CL);

% 
% % Closed Loop 
% % For D controller
% kd = 1;
% OL = tf([0, 10*kd, 0], [1, 5, 5]);
% CL = tf([0, 10*kd, 0], [1, 5+10*kd, 5]);
% figure
% rlocus(OL)
% grid on
% title('Root Locus')
% figure
% step(CL)
% grid on
% xline(1,'--r','Steadx State')
% title('Closed-Loop Step Response')
% p = pole(CL);
% info = stepinfo(CL);



% % Closed Loop 
% % For PI controller
% kp = 200;
% ki = 100; 
% OL = tf([0, 0, 10*kp, 10*ki], [1, 5, 5, 0]);
% CL = tf([0, 0, 10*kp, 10*ki], [1, 5, 5+10*kp, 10*ki]);
% figure
% rlocus(OL)
% grid on
% title('Root Locus')
% figure
% step(CL)
% grid on
% xline(1,'--r','Steadx State')
% title('Closed-Loop Step Response')
% p = pole(CL);
% info = stepinfo(CL);


% % Closed Loop 
% % For PD controller
% kp = 100;
% kd = 100; 
% OL = tf([0, 10*kd, 10*kp], [1, 5, 5]);
% CL = tf([0, 10*kd, 10*kp], [1, 5+10*kd, 5+10*kp]);
% figure
% rlocus(OL)
% grid on
% title('Root Locus')
% figure
% step(CL)
% grid on
% xline(1,'--r','Steadx State')
% title('Closed-Loop Step Response')
% p = pole(CL);
% info = stepinfo(CL);



% Closed Loop 
% For PID controller
kp = 100;
ki = 100;
kd = 100; 
OL = tf([10*kd, 10*kp, 10*ki], [1, 5, 5, 0]);
CL = tf([10*kd, 10*kp, 10*ki], [1, 5+10*kd, 5+10*kp, 10*ki]);
figure
rlocus(OL)
grid on
title('Root Locus')
figure
step(CL)
grid on
xline(1,'--r','Steadx State')
title('Closed-Loop Step Response')
p = pole(CL);
info = stepinfo(CL);










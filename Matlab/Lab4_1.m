clc; 
clear variables; 

vp = [360   0;
      360   120;
      360   -120;
      ];

Zs = 1j*24; Zm = 1j*6;
a1 = 1 * exp(1i * deg2rad(120));
a2 = a1^2;

A = [1  1   1;
     1  a2  a1;
     1  a1  a2;
     ];

vp_mag = vp(:, 1);
vp_angle = vp(:, 2)*(pi/180);

vabc = vp(:, 1) .* exp(1j*vp_angle);

v012 = inv(A)*vabc;

z012 = [Zs+2*Zm    0        0;
        0        Zs - Zm    0;
        0          0        Zs - Zm;
        ];

I012 = inv(z012)*v012;
Iabc = A*I012;

disp("Line to Ground voltages are: ")
disp(vabc)
disp("phase impedances is: ")
disp(Zs)
disp("mutual impedances is: ")
disp(Zm)
disp("The phase currents are: ")
disp(Iabc)
disp("The sequence voltages are: ")
disp(v012)
disp("The sequence currents are: ")
disp(I012)
clc 
clear vars;
%         bus1 bus2  R     X  
line_data = [1  2    0  0.125;
             1  3    0  0.25;
             1  4    0  0.4;
             2  3    0  0.25;
             2  4    0  0.2;
];

nbus = max(max(line_data(:, 1:2)));  % Number of buses
nline = size(line_data, 1);         % Number of Transmission lines

YBus = zeros(nbus, nbus);           % Initialize Y-bus
R = line_data(:, 3);
X = line_data(:, 4);
Z = R + 1j*X;
Y = 1./Z;

for i = 1:nline 
    from = line_data(i,1);
    to = line_data(i,2);
    Yij = Y(i);
    YBus(from, from) = YBus(from, from) + Yij;
    YBus(to, to) =  YBus(to, to) + Yij;
    YBus(from, to) = YBus(from, to) - Yij;
    YBus(to, from) = YBus(to, from) - Yij; 
end 

ZBus = inv(YBus);

disp(' '); 
disp('Calculated Y-bus Admittance Matrix:');
disp(YBus);

disp(' '); 
disp('Calculated Z-bus Admittance Matrix:');
disp(ZBus);




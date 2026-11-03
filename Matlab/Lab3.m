clc 
clear vars;

%bus_data   bus   V            P       Q  
bd =        [1    1 + 1j*0    NaN     NaN; %slack bus 
             2    NaN         -0.2    -0.1; %PQ bus 
             3    NaN         -0.2    -0.1; %PQ bus 
             ]; 

%line_data  bus1 bus2   R     X  
ld =        [1   2      0     0.07;
             1   3      0     0.05;
             2   3      0     0.1
             ];

nbus = max(bd(:, 1));  % Number of buses
nline = size(ld, 1);   % Number of Transmission lines

YBus = zeros(nbus, nbus);           % Initialize Y-bus
R = ld(:, 3);
X = ld(:, 4);
Z = R + 1j*X;
Y = 1./Z;

for i = 1:nline 
    from = ld(i,1);
    to = ld(i,2);
    Yij = Y(i);
    YBus(from, from) = YBus(from, from) + Yij;
    YBus(to, to) =  YBus(to, to) + Yij;
    YBus(from, to) = YBus(from, to) - Yij;
    YBus(to, from) = YBus(to, from) - Yij; 
end 
disp("YBUS: ")
disp(YBus)

%Assuming load bus voltages 
for e = 1:nbus
    if isnan(bd(e, 2))
        bd(e, 2) = 1 + 1j*0;
    end
end

% Number of Iteration 
it = 2; 

%Gauss seidal iteration
for m = 1:it 
    for p = 2:nbus 
        t2 = 0;
            for q = 1:nbus
                if q ~= p
                    t2 = t2 + YBus(p, q)*bd(q, 2);              
                end
            end
        t1 = (bd(p, 3) - 1j*bd(p, 4))/conj(bd(p, 2));
        bd(p, 2) = (1/YBus(p, p))*( t1 - t2);
    end
end

disp("Voltage of each bus is:")
disp(bd(:, 2))




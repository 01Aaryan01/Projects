clc;
clearvars;

f = 50; %  Hz 
V_Level = 345e3; %  V
L = 130e3; %  Meters  
r = 0.036; %  ohm/KM-phase
l = 0.8e-3; %  H/KM-phase 
c = 0.0112e-6; %  F/KM-phase 
S_l = 270e6; %  VA 
phi_l = acos(0.8); 

% Lumped Parameters 
function [R, XL, XC, Z, Y, Zc, g] = lmp_p(r,l,c,f,L)  
    R = r*L; % ohm/phase 
    XL = 2*pi*f*l*L; % ohm/phase 
    XC = 1/(2*pi*f*c*L); % ohm/phase 
    Z = R + 1j*XC; 
    Y = 1j*XC;
    Zc = sqrt(Z/Y);
    g = sqrt(Z*Y); 
end 

[R, XL, XC, Z, Y, Zc, g] = lmp_p(r,l,c,f,L);

% Short TL Model 
function [A,B,C,D] = STL(Z) 
    A = 1; 
    D = 1; 
    B = Z;
    C = 0;
end 

% Medium Nominal pi TL Model 
function [A,B,C,D] = MPi(Z, Y)
    A = 1 + (Y*Z)/2;
    D = A;
    B = Z*(1 + (Y*Z)/4);
    C = Y;
end 

% Medium Nominal T TL Model 
function [A,B,C,D] = MT(Z, Y)
    A = 1 + (Y*Z)/2;
    D = A;
    C = Y*(1 + (Y*Z)/4);
    B = Z;
end 

% Long Rigorous TL Model 
function [A,B,C,D] = LR(Zc, g, L)
    A = cosh(g*L);
    D = A;
    C = (1/Zc)*(sinh(g*L));
    B = Zc*(sinh(g*L));
end

% Long Equivalent pi Method
function [A,B,C,D] = LPi(Z, Y, g, L)
    Z1 = (Z*tanh(g*L))/(g*L);
    Y1 = (Y*tanh((g*L)/2))/((g*L)/2);
    A = 1 + ((Y1*Z1)/2);
    D = A;
    C = Y1*(1 + ((Y1*Z1)/4));
    B = Z1;
end

% Long Equivalent T Method
function [A,B,C,D] = LT(Z, Y, g, L)
    Z1 = (Z*tanh((g*L)/2))/((g*L)/2);
    Y1 = (Y*tanh(g*L))/(g*L);
    A = 1 + ((Y1*Z1)/2);
    D = A;
    C = Y1*(1 + ((Y1*Z1)/4));
    B = Z1;
end


% Model names
model_names = {
    'Short TL Model', 
    'Medium Pi TL Model', 
    'Medium T TL Model', 
    'Long Rigorous TL Model',
    'Long Pi TL Model',
    'Long T TL Model'
};

% Function handles for each model
model_funcs = {
    @() STL(Z),
    @() MPi(Z, Y),
    @() MT(Z, Y),
    @() LR(Zc, g, L),
    @() LPi(Z, Y, g, L),
    @() LT(Z, Y, g, L)
};

% Loop through each model
for i = 1:length(model_names)
    [A, B, C, D] = model_funcs{i}();
    
    fprintf('\n--- %s ---\n', model_names{i});
    fprintf('A = %.4f + %.4fj\n', real(A), imag(A));
    fprintf('B = %.4f + %.4fj\n', real(B), imag(B));
    fprintf('C = %.4f + %.4fj\n', real(C), imag(C));
    fprintf('D = %.4f + %.4fj\n', real(D), imag(D));
end








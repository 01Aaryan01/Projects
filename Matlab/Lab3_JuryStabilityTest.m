%{ 
    # Stabiltiy Analysis of Closed-loop System in Z-plane
    * Jury Stability Test
    - Analysis done based on the coefficients of characteristics equation
      P(z) = a0 * z^n + a1 * z^(n-1) + ... + a(n-1) * z + an
      
      eg: For 4th order P(z)
      P(z) = a0 * z^4 + a1 * z^3 + a2 * z^2 + a3 * z + a4
    
    % Stabitlity test
    1. |an| < a0
    2. P(z)|(z =  1) > 0
    3. For,
       Even n, P(z)|(z = -1) > 0
       Odd  n, P(z)|(z = -1) < 0
    4. From Jury Table
       |b(n-1)| > |b0|
       |b(n-1)| > |b0|
       .
       .
       .
       |q2| > |q0|
%}

% System Parameters
CharEq = [1, -1.2, 0.07, 0.3, -0.08]; % Characteristic Equation in z-domain
% CharEq = [1, -1.3, -0.08, 0.24]; % Characteristic Equation in z-domain
% CharEq = [1, -1.1, -0.1, 0.2]; % Characteristic Equation in z-domain

% Stability Test
%% First Test
% |an| < a0
%{
    eg: For 4th order P(z)
    The characteristics equation:
    P(z) = a0 * z^4 + a1 * z^3 + a2 * z^2 + a3 * z + a4
    So, CharEq = [a0, a1, a2, a3, a4];
%}

an = CharEq(end);
a0 = CharEq(1);
condition_1 = ( abs(an) < abs(a0) ); % Check for 1st condition

%% Second Test
% P(z)|(z =  1) > 0
%{
    P(z) evaluated for z = 1
    eg: For 4th order P(z)
    P(z) = a0 * z^4 + a1 * z^3 + a2 * z^2 + a3 * z + a4
    P(1) = a0 * 1^4 + a1 * 1^3 + a2 * 1^2 + a3 * 1 + a4
    P(1) = a0       + a1       + a2       + a3     + a4
    i.e. simply adding all the coefficient of characteristic equation
%}

Pz1 = sum(CharEq);
condition_2 = (Pz1 > 0); % Check for 2nd condition

%% Third Test
% For,
%   Even n, P(z)|(z = -1) > 0
%   Odd  n, P(z)|(z = -1) < 0

%{
    P(z) evaluated for z = -1
    eg: For 4th order P(z)
    P(z)  = a0 *    z^4 + a1 *    z^3 + a2 *    z^2 + a3 *    z + a4
    P(-1) = a0 * (-1)^4 + a1 * (-1)^3 + a2 * (-1)^2 + a3 * (-1) + a4
    P(-1) = a0          - a1          + a2          - a3        + a4

    Implementation in matlab
    P(-1) = [a0, a1, a2, a3, a4] .* [     1,     -1,      1,     -1,      1]
    P(-1) = [a0, a1, a2, a3, a4] .* [(-1)^4, (-1)^3, (-1)^2, (-1)^1, (-1)^0]
    P(-1) = [a0, a1, a2, a3, a4] .* (-1).^[4,3,2,1,0]
%}

order = length(CharEq) - 1; % Power of z is from n to 0
vector_of_minus1 = (-1).^(order:-1:0); % Element-wise power operation
Pzm1 = sum(CharEq .* vector_of_minus1);

if rem(order,2) == 0 % Check for even n
    condition_3 = (Pzm1 > 0);
else             % For Odd n
    condition_3 = (Pzm1 < 0);
end

%% Fourth Test
%{
    From Jury Table
    |b(n-1)| > |b0|
    |b(n-1)| > |b0|
    .
    .
    .
    |q2| > |q0|

    eg: For 4th order P(z)
    P(z) = a0 * z^4 + a1 * z^3 + a2 * z^2 + a3 * z + a4
    So, CharEq = [a0, a1, a2, a3, a4];
    First  row: a4, a3, a2, a1, a0
    Second row: a0, a1, a2, a3, a4
%}

% Jury Table
condition_4 = true; % Initilize with stable condition

if order > 2 % Jury Table is only required for system higher than 2nd order
    JuryTable = NaN(2*order-3,(order+1)); % Matrix with NaN value of size (2*n-3) by (n+1)
    
    % Insert characteristic equation in the Jury Table
    JuryTable(1,:) = flip(CharEq); % Flip's the vector CharEq
    JuryTable(2,:) = CharEq;
    
    % Initilize indexing variables
    row_index = 3; % Row 1 and row 2 are aldreay defined, we start from row 3
    column_count = order+1; % (n+1) to account for z^0 coefficient 
    
    for j = row_index:2:(2*order-3) % Going through all rows with step of 2
        
        temp_matrix = JuryTable((j-2:j-1),1:column_count); % Temporay matrix to store value of previous 2 rows
        first_column = temp_matrix(:,1);
        
        column_count = column_count - 1; % Column for new row is one less than previous
        temp_vector = NaN(1, column_count);
        
        for k = 1:column_count % For going through all column
            kth_column = temp_matrix(:,k+1);
            temp_vector(k) = det([first_column, kth_column]);
        end
        
        % Check 4th condition
        last_element = temp_vector(end);
        first_element = temp_vector(1);
        temp_condition = abs(last_element) > abs(first_element);
        condition_4 = condition_4 && temp_condition;
        
        if column_count == 3
            JuryTable(end,1:column_count) = flip(temp_vector);
        else
            JuryTable(  j, 1:column_count) = flip(temp_vector);
            JuryTable(j+1, 1:column_count) = temp_vector;
        end
    end
end

%% Print stability result

% If all condition are true, then the system is stable
if condition_1 && condition_2 && condition_3 && condition_4
    fprintf('System is stable.\n');
else
    fprintf('System is unstable.\n');
end
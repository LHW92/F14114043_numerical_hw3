%% Q1
clc; clear; close all;

% 已知數據
x = [0.698 0.733 0.768 0.803]; 
y = [0.7661 0.7432 0.7193 0.6946];

x_tar = 0.75; % 目標插值點

fprintf('Degree 1 Approximation:\n');
fprintf('i, j\t L0\t\t\t L1\t\t\t y_tar\n');
fprintf('-------------------------------------------------\n');

n = length(x); % 數據點數量

for i = 1:n-1
    for j = i+1:n
        % 計算 Lagrange 基函數 L0, L1
        L0 = (x_tar - x(j)) / (x(i) - x(j));
        L1 = (x_tar - x(i)) / (x(j) - x(i));
        
        % 計算插值結果 y_tar
        y_tar = L0 * y(i) + L1 * y(j);
        
        % 輸出結果
        fprintf('%d, %d\t %.6f\t %.6f\t %.6f\n', i-1, j-1, L0, L1, y_tar);
    end
end

fprintf('\nDegree 2 Approximation:\n');
fprintf('i, j, k\t L0\t\t\t L1\t\t\t L2\t\t\t y_tar\n');
fprintf('------------------------------------------------------------\n');

for i = 1:n-2
    for j = i+1:n-1
        for k = j+1:n
            % 計算 Lagrange 基函數 L0, L1, L2
            L0 = ((x_tar - x(j)) * (x_tar - x(k))) / ((x(i) - x(j)) * (x(i) - x(k)));
            L1 = ((x_tar - x(i)) * (x_tar - x(k))) / ((x(j) - x(i)) * (x(j) - x(k)));
            L2 = ((x_tar - x(i)) * (x_tar - x(j))) / ((x(k) - x(i)) * (x(k) - x(j)));

            % 計算插值結果 y_tar
            y_tar = L0 * y(i) + L1 * y(j) + L2 * y(k);

            % 輸出結果
            fprintf('%d, %d, %d\t %.6f\t %.6f\t %.6f\t %.6f\n', i-1, j-1, k-1, L0, L1, L2, y_tar);
        end
    end
end

fprintf('\nDegree 3 Approximation:\n');
fprintf('i, j, k, m\t L0\t\t\t L1\t\t\t L2\t\t\t L3\t\t\t y_tar\n');
fprintf('------------------------------------------------------------------------\n');

for i = 1:n-3
    for j = i+1:n-2
        for k = j+1:n-1
            for m = k+1:n
                % 計算 Lagrange 基函數 L0, L1, L2, L3
                L0 = ((x_tar - x(j)) * (x_tar - x(k)) * (x_tar - x(m))) / ((x(i) - x(j)) * (x(i) - x(k)) * (x(i) - x(m)));
                L1 = ((x_tar - x(i)) * (x_tar - x(k)) * (x_tar - x(m))) / ((x(j) - x(i)) * (x(j) - x(k)) * (x(j) - x(m)));
                L2 = ((x_tar - x(i)) * (x_tar - x(j)) * (x_tar - x(m))) / ((x(k) - x(i)) * (x(k) - x(j)) * (x(k) - x(m)));
                L3 = ((x_tar - x(i)) * (x_tar - x(j)) * (x_tar - x(k))) / ((x(m) - x(i)) * (x(m) - x(j)) * (x(m) - x(k)));

                % 計算插值結果 y_tar
                y_tar = L0 * y(i) + L1 * y(j) + L2 * y(k) + L3 * y(m);

                % 輸出結果
                fprintf('%d, %d, %d, %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', i-1, j-1, k-1, m-1, L0, L1, L2, L3, y_tar);
            end
        end
    end
end

%% Q2
clc; clear; close all;

x0=0.3; x1=0.4; x2=0.5; x3=0.6;
y0=-0.440818; y1=-0.270320; y2=-0.106531; y3=0.051188;  % y = x - e^(-x)

y=0;

x=(((y-y1)*(y-y2)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3)))*x0+(((y-y0)*(y-y2)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3)))*x1+(((y-y0)*(y-y1)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3)))*x2+(((y-y0)*(y-y1)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2)))*x3;
y=exp(-x);

fprintf('Approximate root: x ≈ %.6f\n', x);

%% Q3 (use method in textbook)

clc; clear;

% Given data (all 5 points)
T = [0, 3, 5, 8, 13];  % Time data points
D = [0, 200, 375, 620, 990]; % Distance (position)
V = [75, 77, 80, 74, 72]; % Speed (velocity)

% Given data (last 4 points)
% T = [3, 5, 8, 13];  % Time data points
% D = [200, 375, 620, 990]; % Distance (position)
% V = [77, 80, 74, 72]; % Speed (velocity)

% Given data (last 3 points)
% T = [5, 8, 13];  % Time data points
% D = [375, 620, 990]; % Distance (position)
% V = [80, 74, 72]; % Speed (velocity)

t_target = 10; % The time at which we evaluate the basis functions


L_values = zeros(1, length(T));
L_deriv_values = zeros(1, length(T));
for i = 1:length(T)
    L_values(i) = lagrange_basis(T, i, t_target);              % Compute Lagrange basis functions
    L_deriv_values(i) = lagrange_basis_derivative(T, i, T(i)); % Compute derivatives of Lagrange basis functions
end

% Compute interpolated position
predicted_position = hermite_interpolation(T, D, V, t_target); 

% Compute interpolated speed
predicted_speed = hermite_interpolation(T, V, zeros(1, length(V)), t_target);

fprintf("(a) Predicted Position at t = %.2f s: %.2f ft\n", t_target, predicted_position);
fprintf("    Predicted Speed at t = %.2f s: %.2f ft/s\n", t_target, predicted_speed);

% Display results
% fprintf('Lagrange basis function values at t = %.2f:\n', t_target);
% for i = 1:length(T)
%     fprintf('L_%d(%.2f) = %.6f\n', i-1, t_target, L_values(i));
% end

% Display results
%fprintf('Lagrange basis function derivatives at t = %.2f:\n', t_target);
% for i = 1:length(T)
%     fprintf("L_%d'(%.2f) = %.6f\n", i-1, T(i), L_deriv_values(i));
% end

% % Compute Hermite basis functions for each j
% H_values = zeros(1, length(T));
% for j = 1:length(T)
%     H_values(j) = hermite_basis_function(T, j, t_target);
% end
% 
% % Display results
% fprintf('Hermite basis functions at t = %.2f:\n', t_target);
% for j = 1:length(T)
%     fprintf("H_%d(%.2f) = %.6f\n", j-1, t_target, H_values(j));
% end
% 
% % Compute Hermite auxiliary basis functions for each j
% Hhat_values = zeros(1, length(T));
% for j = 1:length(T)
%     Hhat_values(j) = hermite_auxiliary_basis(T, j, t_target);
% end
% 
% % Display results
% fprintf('Hermite auxiliary basis functions at t = %.2f:\n', t_target);
% for j = 1:length(T)
%     fprintf("Hhat_%d(%.2f) = %.6f\n", j-1, t_target, Hhat_values(j));
% end

function L = lagrange_basis(T, i, t)
    % Compute the Lagrange basis function L_i(t)
    n = length(T); % Number of data points
    L = 1;
    for j = 1:n
        if j ~= i
            L = L * (t - T(j)) / (T(i) - T(j));
        end
    end
end

function L_deriv = lagrange_basis_derivative(T, i, t)
    % Compute the derivative of the Lagrange basis function L_i'(t)
    n = length(T); % Number of data points
    L_deriv = 0;
    
    for k = 1:n
        if k ~= i
            term = 1 / (T(i) - T(k));
            prod_term = 1;
            for j = 1:n
                if j ~= i && j ~= k
                    prod_term = prod_term * (t - T(j)) / (T(i) - T(j));
                end
            end
            L_deriv = L_deriv + term * prod_term;
        end
    end
end

% function H = hermite_basis_function(T, j, x)
%     % Compute the Hermite basis function H_{n,j}(x)
% 
%     % Compute Lagrange basis polynomial L_{n,j}(x)
%     L = lagrange_basis(T, j, x);
% 
%     % Compute derivative L_{n,j}'(x_j)
%     L_prime = lagrange_basis_derivative(T, j, T(j));
% 
%     % Compute Hermite basis function
%     H = (1 - 2 * (x - T(j)) * L_prime) * (L^2);
% end
% 
% function Hhat = hermite_auxiliary_basis(T, j, x)
%     % Compute the auxiliary Hermite basis function Hhat_{n,j}(x)
% 
%     % Compute Lagrange basis polynomial L_{n,j}(x)
%     L = lagrange_basis(T, j, x);
% 
%     % Compute auxiliary Hermite basis function
%     Hhat = (x - T(j)) * (L^2);
% end

function H = hermite_interpolation(T, F, F_deriv, x)
    % Compute the Hermite interpolating polynomial H_{2n+1}(x)
    n = length(T);
    H = 0;
    
    for j = 1:n
        L = lagrange_basis(T, j, x);
        L_deriv = lagrange_basis_derivative(T, j, T(j));

        H_j = (1 - 2*(x - T(j)) * L_deriv) * (L^2);
        Hhat_j = (x - T(j)) * (L^2);
        
        H = H + F(j) * H_j + F_deriv(j) * Hhat_j;
    end
end

t_range = linspace(0, 13, 1000); % Fine time grid

% Compute velocity using Hermite interpolation derivative
V_predicted = arrayfun(@(t) hermite_interpolation_derivative(T, D, V, t), t_range);

% Find the first time when V(t) > 80.67 ft/s
speed_limit = 80.67;
exceed_idx = find(V_predicted > speed_limit, 1);
if isempty(exceed_idx)
    fprintf("(b) The car never exceeds 55 mi/h.\n");
else
    first_exceed_time = t_range(exceed_idx);
    fprintf("(b) The car first exceeds 55 mi/h at t = %.2f s\n", first_exceed_time);
end

% Find the maximum predicted speed
[max_speed, max_idx] = max(V_predicted);
max_speed_time = t_range(max_idx);
fprintf("(c) The predicted maximum speed is %.2f ft/s at t = %.2f s\n", max_speed, max_speed_time);

% Function for Hermite interpolation derivative
function H_deriv = hermite_interpolation_derivative(T, F, F_deriv, x)
    % Compute the derivative of the Hermite interpolating polynomial H'(x)
    n = length(T);
    H_deriv = 0;
    
    for j = 1:n
        L = lagrange_basis(T, j, x);
        L_deriv = lagrange_basis_derivative(T, j, x);
        L_prime_at_Tj = lagrange_basis_derivative(T, j, T(j));

        H_j_deriv = 2 * L * L_deriv * (1 - 2 * (x - T(j)) * L_prime_at_Tj) - 2 * L_prime_at_Tj * L^2;
        Hhat_j_deriv = L^2 + 2 * (x - T(j)) * L * L_deriv;
        
        H_deriv = H_deriv + F(j) * H_j_deriv + F_deriv(j) * Hhat_j_deriv;
    end
end

%% Q3 (use pchip)

clc; clear;

% Given Data
T = [0, 3, 5, 8, 13];   % Time (seconds)
D = [0, 200, 375, 620, 990]; % Distance (feet)
V = [75, 77, 80, 74, 72]; % Speed (feet/sec)

% Target time
t_target = 10;

% Hermite Interpolation using PCHIP (Piecewise Cubic Hermite Interpolating Polynomial)
D_interp = @(t) pchip(T, D, t);
V_interp = @(t) pchip(T, V, t); % Speed is the derivative

% (a) Predict position and speed at t = 10
D_10 = D_interp(t_target);
V_10 = V_interp(t_target);

fprintf('(a) Predicted position at t = 10s: %.2f feet\n', D_10);
fprintf('    Predicted speed at t = 10s: %.2f feet/sec\n', V_10);

% (b) Check if car exceeds 55 mi/h speed limit
mph_to_ftps = 22/15; % 1 mph = 1.467 ft/sec
speed_limit_ftps = 55 * mph_to_ftps;

t_fine = linspace(0, 13, 100); % More points for analysis
V_fine = V_interp(t_fine);

exceeds_limit = find(V_fine > speed_limit_ftps, 1, 'first');

if isempty(exceeds_limit)
    fprintf('(b) Car never exceeds 55 mph.\n');
else
    fprintf('(b) Car first exceeds 55 mph at t = %.2f seconds.\n', t_fine(exceeds_limit));
end

% (c) Predicted maximum speed
[max_speed, max_idx] = max(V_fine);
fprintf('(c) Predicted maximum speed: %.2f feet/sec at t = %.2f seconds\n', max_speed, t_fine(max_idx));
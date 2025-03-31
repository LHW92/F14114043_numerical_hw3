%% use method in textbook
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
%% use pchip

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
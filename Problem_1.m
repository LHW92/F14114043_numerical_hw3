clc; clear; close all;

% 已知數據
x = [0.698 0.733 0.768 0.803]; 
y = [0.7661 0.7432 0.7193 0.6946];

x_tar = 0.75; % 目標插值點
true_value = cos(x_tar); % 真實值 cos(0.75) = 0.7317

fprintf('Degree 1 Approximation:\n');
fprintf('i, j\t L0\t\t\t L1\t\t\t y_tar\t\t Error Bound\n');
fprintf('-------------------------------------------------------------\n');

n = length(x); % 數據點數量

for i = 1:n-1
    for j = i+1:n
        % 計算 Lagrange 基函數 L0, L1
        L0 = (x_tar - x(j)) / (x(i) - x(j));
        L1 = (x_tar - x(i)) / (x(j) - x(i));
        
        % 計算插值結果 y_tar
        y_tar = L0 * y(i) + L1 * y(j);
        
        % 計算誤差界限
        error_bound = (1 / 2) * abs((x_tar - x(i)) * (x_tar - x(j)));
        
        % 輸出結果
        fprintf('%d, %d\t %.6f\t %.6f\t %.6f\t %.6f\n', i-1, j-1, L0, L1, y_tar, error_bound);
    end
end

fprintf('\nDegree 2 Approximation:\n');
fprintf('i, j, k\t L0\t\t\t L1\t\t\t L2\t\t\t y_tar\t\t Error Bound\n');
fprintf('------------------------------------------------------------------\n');

for i = 1:n-2
    for j = i+1:n-1
        for k = j+1:n
            % 計算 Lagrange 基函數
            L0 = ((x_tar - x(j)) * (x_tar - x(k))) / ((x(i) - x(j)) * (x(i) - x(k)));
            L1 = ((x_tar - x(i)) * (x_tar - x(k))) / ((x(j) - x(i)) * (x(j) - x(k)));
            L2 = ((x_tar - x(i)) * (x_tar - x(j))) / ((x(k) - x(i)) * (x(k) - x(j)));
            
            % 計算插值結果 y_tar
            y_tar = L0 * y(i) + L1 * y(j) + L2 * y(k);
            
            % 計算誤差界限
            error_bound = (1 / 6) * abs((x_tar - x(i)) * (x_tar - x(j)) * (x_tar - x(k)));
            
            % 輸出結果
            fprintf('%d, %d, %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', i-1, j-1, k-1, L0, L1, L2, y_tar, error_bound);
        end
    end
end

fprintf('\nDegree 3 Approximation:\n');
fprintf('i, j, k, m\t L0\t\t\t L1\t\t\t L2\t\t\t L3\t\t\t y_tar\t\t Error Bound\n');
fprintf('------------------------------------------------------------------------------------\n');

for i = 1:n-3
    for j = i+1:n-2
        for k = j+1:n-1
            for m = k+1:n
                % 計算 Lagrange 基函數
                L0 = ((x_tar - x(j)) * (x_tar - x(k)) * (x_tar - x(m))) / ((x(i) - x(j)) * (x(i) - x(k)) * (x(i) - x(m)));
                L1 = ((x_tar - x(i)) * (x_tar - x(k)) * (x_tar - x(m))) / ((x(j) - x(i)) * (x(j) - x(k)) * (x(j) - x(m)));
                L2 = ((x_tar - x(i)) * (x_tar - x(j)) * (x_tar - x(m))) / ((x(k) - x(i)) * (x(k) - x(j)) * (x(k) - x(m)));
                L3 = ((x_tar - x(i)) * (x_tar - x(j)) * (x_tar - x(k))) / ((x(m) - x(i)) * (x(m) - x(j)) * (x(m) - x(k)));
                
                % 計算插值結果 y_tar
                y_tar = L0 * y(i) + L1 * y(j) + L2 * y(k) + L3 * y(m);
                
                % 計算誤差界限
                error_bound = (1 / 24) * abs((x_tar - x(i)) * (x_tar - x(j)) * (x_tar - x(k)) * (x_tar - x(m)));
                
                % 輸出結果
                fprintf('%d, %d, %d, %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', i-1, j-1, k-1, m-1, L0, L1, L2, L3, y_tar, error_bound);
            end
        end
    end
end

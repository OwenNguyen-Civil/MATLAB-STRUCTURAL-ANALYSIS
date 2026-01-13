clc; clear; close all;
fprintf('=======================================================\n');
fprintf('      GENERAL STRESS CALCULATION PROGRAM (V4.2)        \n');
fprintf('     * User Input Mode - Unsymmetrical Bending * \n');
fprintf('=======================================================\n');

% =========================================================
% 1. NHẬP TẢI TRỌNG (INPUT LOADS)
% =========================================================
fprintf('\n--- BƯỚC 1: NHẬP TẢI TRỌNG ---\n');
Nz = input('1. Nhập Lực dọc Nz (N) [Dương = Kéo, Âm = Nén]: ');
Mx_Nm = input('2. Nhập Mô-men Mx (N.m) [Dương = Căng thớ dưới]: ');
My_Nm = input('3. Nhập Mô-men My (N.m) [Uốn quanh trục đứng]: ');

Mx = Mx_Nm * 1000; % Đổi sang N.mm
My = My_Nm * 1000; % Đổi sang N.mm

% =========================================================
% 2. NHẬP TỌA ĐỘ ĐỈNH TIẾT DIỆN (INPUT VERTEX COORDINATES)
% =========================================================
fprintf('\n--- BƯỚC 2: NHẬP TỌA ĐỘ CÁC ĐỈNH [x y] ---\n');
fprintf('Ví dụ nhập cho thép L: [0 0; 100 0; 100 10; 10 10; 10 100; 0 100]\n');
poly_matrix = input('>> Nhập ma trận tọa độ: ');

try
    x_poly = poly_matrix(:, 1)'; 
    y_poly = poly_matrix(:, 2)'; 
catch
    error('LỖI ĐỊNH DẠNG: Vui lòng nhập ma trận 2 cột [x y]!');
end

% =========================================================
% 3. NHẬP ĐIỂM CẦN TÍNH ỨNG SUẤT (INPUT SPECIFIC POINTS)
% =========================================================
fprintf('\n--- BƯỚC 3: NHẬP CÁC ĐIỂM CẦN TÍNH [x y] ---\n');
pts_matrix = input('>> Nhập ma trận điểm tính toán: ');

try
    x_pts = pts_matrix(:, 1)';
    y_pts = pts_matrix(:, 2)';
catch
    error('LỖI ĐỊNH DẠNG: Vui lòng nhập ma trận 2 cột [x y]!');
end

% =========================================================
% 4. XỬ LÝ HÌNH HỌC (GEOMETRY PROCESSING)
% =========================================================
pgon = polyshape(x_poly, y_poly);
Area = area(pgon);          
[Cx, Cy] = centroid(pgon);  

Ixx_O = 0; Iyy_O = 0; Ixy_O = 0; 
n = length(x_poly);

% Khép kín vòng lặp tọa độ để tính tích phân đường
if x_poly(1) ~= x_poly(end) || y_poly(1) ~= y_poly(end)
    x_loop = [x_poly, x_poly(1)]; y_loop = [y_poly, y_poly(1)]; n_edge = n;
else
    x_loop = x_poly; y_loop = y_poly; n_edge = n - 1;
end

% Định lý Green (Shoelace formula) tính quán tính
for i = 1:n_edge
    Xi = x_loop(i); Xi1 = x_loop(i+1); Yi = y_loop(i); Yi1 = y_loop(i+1);
    D = Xi*Yi1 - Xi1*Yi;
    
    Ixx_O = Ixx_O + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
    Iyy_O = Iyy_O + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
    Ixy_O = Ixy_O + (Xi*Yi1 + 2*Xi*Yi + 2*Xi1*Yi1 + Xi1*Yi)*D/24;
end

% Dời về trọng tâm (Parallel Axis Theorem)
Ix = Ixx_O - Area * Cy^2; 
Iy = Iyy_O - Area * Cx^2; 
Ixy = Ixy_O - Area * Cx * Cy; 

% --- Tính góc xoay trục chính (Dùng atan2 để tránh lỗi chia cho 0) ---
numerator = -2 * Ixy;
denominator = Ix - Iy;
alpha_rad = 0.5 * atan2(numerator, denominator); 
alpha_deg = rad2deg(alpha_rad);

% =========================================================
% 5. TÍNH TOÁN & BÁO CÁO (CALCULATION & REPORT)
% =========================================================
fprintf('\n=================================================\n');
fprintf('             KẾT QUẢ TÍNH TOÁN (V4.2)            \n');
fprintf('=================================================\n');
fprintf('1. Đặc trưng hình học:\n');
fprintf('   - Diện tích A:        %.2f mm^2\n', Area);
fprintf('   - Trọng tâm C:        (%.2f, %.2f) mm\n', Cx, Cy);
fprintf('   - Ix (Ngang):         %.3e mm^4\n', Ix);
fprintf('   - Iy (Đứng):          %.3e mm^4\n', Iy);
fprintf('   - Ixy (Ly tâm):       %.3e mm^4\n', Ixy);
if abs(Ixy) > 1e-3
    fprintf('      -> Tiết diện KHÔNG ĐỐI XỨNG (Uốn xiên).\n');
    fprintf('      -> Góc xoay trục chính: %.2f độ\n', alpha_deg);
else
    fprintf('      -> Tiết diện ĐỐI XỨNG qua trục XY.\n');
end
fprintf('-------------------------------------------------\n');
fprintf('2. Ứng suất tại các điểm chỉ định (MPa):\n');

Points_X_new = x_pts - Cx;
Points_Y_new = y_pts - Cy;
Stress_Specific = zeros(1, length(x_pts)); 

Det_I = Ix * Iy - Ixy^2; 

for k = 1:length(x_pts)
    x = Points_X_new(k);
    y = Points_Y_new(k);
    
    % [CÔNG THỨC TỔNG QUÁT]
    % Mx dương gây NÉN thớ trên (+y) -> thêm dấu trừ cho thành phần y
    term_N = Nz / Area;
    term_Mx = -(Mx * Iy - My * Ixy) * y / Det_I; 
    term_My = (My * Ix - Mx * Ixy) * x / Det_I;
    
    val = term_N + term_Mx + term_My;
    
    Stress_Specific(k) = val;
    fprintf('Điểm %d (X=%.0f, Y=%.0f) | Sigma = %8.3f MPa\n', k, x_pts(k), y_pts(k), val);
end

% =========================================================
% 6. VẼ BIỂU ĐỒ NHIỆT & ĐƯỜNG TRUNG HÒA
% =========================================================
figure('Name', 'Stress Map V4.2', 'Color', 'k');

% Tạo lưới điểm để vẽ màu
pgon_shifted = translate(pgon, -Cx, -Cy);
[xlim_b, ylim_b] = boundingbox(pgon_shifted);
res = 200; 
xg = linspace(xlim_b(1), xlim_b(2), res);
yg = linspace(ylim_b(1), ylim_b(2), res);
[X_grid, Y_grid] = meshgrid(xg, yg);
in = isinterior(pgon_shifted, X_grid(:), Y_grid(:));
in = reshape(in, size(X_grid));

Stress_Grid = nan(size(X_grid)); 

% Tính ứng suất cho từng điểm trên lưới
Stress_Grid(in) = (Nz / Area) + ...
                  (-(Mx * Iy - My * Ixy) .* Y_grid(in) ./ Det_I) + ...
                  ((My * Ix - Mx * Ixy) .* X_grid(in) ./ Det_I);

% Vẽ Contour
contourf(X_grid, Y_grid, Stress_Grid, 50, 'LineStyle', 'none'); 
hold on; 
plot(pgon_shifted, 'FaceColor', 'none', 'EdgeColor', 'w', 'LineWidth', 2);

% --- VẼ ĐƯỜNG TRUNG HÒA (Neutral Axis - NA) ---
% Phương trình đường trung hòa: Sigma = 0
num_NA = (My * Ix - Mx * Ixy);
den_NA = (Mx * Iy - My * Ixy);

if abs(den_NA) > 1e-10
    slope_NA = num_NA / den_NA;
    x_line = linspace(xlim_b(1), xlim_b(2), 10);
    y_line = slope_NA * x_line;
    plot(x_line, y_line, 'm--', 'LineWidth', 2, 'DisplayName', 'Neutral Axis');
else
    xline(0, 'm--', 'LineWidth', 2, 'DisplayName', 'Neutral Axis');
end

% Vẽ trục hình học
xline(0, 'w-.', 'Alpha', 0.5); yline(0, 'w-.', 'Alpha', 0.5);

% Đánh dấu điểm tính toán
scatter(Points_X_new, Points_Y_new, 80, 'r', 'filled', 'MarkerEdgeColor', 'w');

% Trang trí
colormap(jet); 
c = colorbar; c.Label.String = 'Stress (MPa)'; c.Color = 'w'; c.Label.Color = 'w';
title(sprintf('STRESS DISTRIBUTION\nNeutral Axis Rotation: %.1f deg', rad2deg(atan(num_NA/den_NA))), 'Color', 'w');
axis equal; grid on;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');

% Ghi chú giá trị lên hình
for k = 1:length(x_pts)
    text(Points_X_new(k), Points_Y_new(k), sprintf('  %.1f MPa', Stress_Specific(k)), ...
         'Color', 'y', 'FontSize', 11, 'FontWeight', 'bold');
end
legend('Location','best','TextColor','w','Color','none');
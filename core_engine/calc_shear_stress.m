function [Tau_max, Debug_Data] = calc_shear_stress(PolyCoords, Vy, Geom_Props)
% CALC_SHEAR_STRESS: Tính ứng suất tiếp bằng phương pháp "Máy in" (Pixel)
% -----------------------------------------------------------------------
% INPUTS:
%   PolyCoords : Ma trận [x y] đỉnh tiết diện
%   Vy         : Lực cắt phương đứng (N)
%   Geom_Props : Struct chứa đặc trưng hình học (lấy từ hàm trước)
%                (Cần: Cx, Cy, Ix)
%
% OUTPUTS:
%   Tau_max    : Ứng suất tiếp lớn nhất (MPa)
%   Debug_Data : Struct chứa dữ liệu để vẽ biểu đồ (Y_grid, Tau_profile)
% -----------------------------------------------------------------------

    % 1. CẤU HÌNH "MÁY IN" (RESOLUTION)
    res_y = 500; % Độ phân giải 500 lát cắt
    y_poly = PolyCoords(:, 2);
    x_poly = PolyCoords(:, 1);
    
    min_y = min(y_poly); 
    max_y = max(y_poly);
    
    % Tạo các lát cắt y (từ dưới lên trên)
    y_slices = linspace(min_y, max_y, res_y)';
    dy = y_slices(2) - y_slices(1);
    
    % Tạo lưới ngang (để đếm bề rộng b)
    % Lấy rộng hơn biên một chút để bao trọn hình
    min_x = min(x_poly); max_x = max(x_poly);
    x_grid_line = linspace(min_x, max_x, 200); 
    dx = x_grid_line(2) - x_grid_line(1);
    
    [X_grid, Y_grid] = meshgrid(x_grid_line, y_slices);
    
    % 2. QUÉT HÌNH (PIXELIZATION)
    pgon = polyshape(x_poly, y_poly);
    in = isinterior(pgon, X_grid(:), Y_grid(:));
    
    % Chuyển vector 'in' thành ma trận đúng kích thước lưới
    Mask = reshape(in, size(X_grid));
    
    % 3. TÍNH TOÁN CƠ HỌC (PHYSICS LOOP)
    % a. Tính bề rộng cắt b(y) tại từng độ cao
    % b(y) = Số lượng pixel ngang * dx
    b_profile = sum(Mask, 2) * dx;
    
    % b. Tính Momen tĩnh Sx(y)
    % Sx = Tích phân (y - Cy) * dA từ y đến đỉnh
    % dA của một hàng = b_profile * dy
    
    Y_dist = y_slices - Geom_Props.Cy; % Khoảng cách đến trục trung hòa
    Area_row = b_profile * dy;         % Diện tích từng lát mỏng
    
    dSx = Y_dist .* Area_row;          % Momen tĩnh của từng lát
    
    % Cộng dồn từ trên xuống (Flip -> Cumsum -> Flip)
    Sx_profile = flipud(cumsum(flipud(dSx)));
    
    % 4. TÍNH ỨNG SUẤT TIẾP (FORMULA: Tau = V*Sx / I*b)
    Tau_profile = zeros(size(Sx_profile));
    
    % Chỉ tính tại những chỗ có vật liệu (b > 0) để tránh chia cho 0
    valid_idx = b_profile > 1e-6; 
    
    % Công thức Zhuravskii
    Tau_profile(valid_idx) = (abs(Vy) * abs(Sx_profile(valid_idx))) ./ ...
                             (Geom_Props.Ix * b_profile(valid_idx));
                             
    % 5. KẾT QUẢ
    Tau_max = max(Tau_profile);
    
    % Đóng gói dữ liệu để vẽ (nếu cần)
    Debug_Data.Y = y_slices;
    Debug_Data.Sx = Sx_profile;
    Debug_Data.b = b_profile;
    Debug_Data.Tau = Tau_profile;
end
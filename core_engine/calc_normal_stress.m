function [Stress_Res, Geom_Props] = calc_normal_stress(PolyCoords, CheckPoints, Loads)
% CALC_NORMAL_STRESS (V3 - FINAL)
% - Tự động quét đỉnh tìm Max/Min
% - Tính góc nghiêng trục trung hòa (Neutral Axis Rotation)
% -----------------------------------------------------------------------
    % Trước khi tạo polyshape, lọc các điểm trùng nhau
    PolyCoords = unique(PolyCoords, 'rows', 'stable');
    % --- 1. TÍNH ĐẶC TRƯNG HÌNH HỌC ---
    x_poly = PolyCoords(:, 1)';
    y_poly = PolyCoords(:, 2)';
    
    if x_poly(1) ~= x_poly(end) || y_poly(1) ~= y_poly(end)
        x_poly = [x_poly, x_poly(1)];
        y_poly = [y_poly, y_poly(1)];
    end
    
    n = length(x_poly) - 1;
    Area = 0; Cx_temp = 0; Cy_temp = 0;
    Ixx_O = 0; Iyy_O = 0; Ixy_O = 0;
    
    for i = 1:n
        Xi = x_poly(i); Xi1 = x_poly(i+1);
        Yi = y_poly(i); Yi1 = y_poly(i+1);
        D = Xi*Yi1 - Xi1*Yi;
        Area = Area + D/2;
        Cx_temp = Cx_temp + (Xi + Xi1)*D/6;
        Cy_temp = Cy_temp + (Yi + Yi1)*D/6;
        Ixx_O = Ixx_O + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
        Iyy_O = Iyy_O + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
        Ixy_O = Ixy_O + (Xi*Yi1 + 2*Xi*Yi + 2*Xi1*Yi1 + Xi1*Yi)*D/24;
    end
    
    Cx = Cx_temp / Area;
    Cy = Cy_temp / Area;
    Ix = Ixx_O - Area * Cy^2;
    Iy = Iyy_O - Area * Cx^2;
    Ixy = Ixy_O - Area * Cx * Cy;
    
    % Góc xoay trục quán tính chính (Alpha) - Chỉ phụ thuộc hình học
    alpha_rad = 0.5 * atan2(-2 * Ixy, Ix - Iy);
    
    % --- 2. TÍNH TOÁN LIÊN QUAN ĐẾN TẢI TRỌNG & TRỤC TRUNG HÒA ---
    Nz = Loads.Nz; Mx = Loads.Mx; My = Loads.My;
    Det_I = Ix * Iy - Ixy^2;
    
    % Hệ số trong phương trình ứng suất: Sigma = Term_N + K1*y + K2*x
    % K1 là hệ số của y, K2 là hệ số của x
    K1 = -(Mx * Iy - My * Ixy) / Det_I; 
    K2 =  (My * Ix - Mx * Ixy) / Det_I;
    Term_N = Nz / Area;
    
    calc_sigma = @(x, y) Term_N + K1 * y + K2 * x;

    % --- TÍNH GÓC TRỤC TRUNG HÒA (BETA) ---
    % Phương trình trục trung hòa: Sigma = 0 => Term_N + K1*y + K2*x = 0
    % => y = (-K2/K1)*x - (Term_N/K1)
    % Hệ số góc (Slope) = -K2/K1
    
    if abs(K1) < 1e-10 % Tránh chia cho 0 (khi uốn quanh trục yếu nhất)
        Beta_deg = 90; % Trục thẳng đứng
    else
        Slope_NA = -K2 / K1;
        Beta_rad = atan(Slope_NA);
        Beta_deg = rad2deg(Beta_rad);
    end

    % Đóng gói kết quả hình học
    Geom_Props = struct('Area', Area, 'Cx', Cx, 'Cy', Cy, ...
                        'Ix', Ix, 'Iy', Iy, 'Ixy', Ixy, ...
                        'Alpha_Principal_deg', rad2deg(alpha_rad), ...
                        'Beta_NeutralAxis_deg', Beta_deg);

    % --- 3. TÍNH TẠI CHECKPOINTS ---
    if ~isempty(CheckPoints)
        x_chk = CheckPoints(:, 1)' - Cx;
        y_chk = CheckPoints(:, 2)' - Cy;
        Stress_Res = zeros(length(x_chk), 1);
        for k = 1:length(x_chk)
            Stress_Res(k) = calc_sigma(x_chk(k), y_chk(k));
        end
    else
        Stress_Res = [];
    end
    
    % --- 4. AUTO SCAN MAX/MIN (GLOBAL) ---
    x_verts = PolyCoords(:, 1) - Cx;
    y_verts = PolyCoords(:, 2) - Cy;
    All_Vert_Stress = zeros(length(x_verts), 1);
    for m = 1:length(x_verts)
        All_Vert_Stress(m) = calc_sigma(x_verts(m), y_verts(m));
    end
    
    Real_Max = max(All_Vert_Stress);
    Real_Min = min(All_Vert_Stress);

    % --- 5. BÁO CÁO KẾT QUẢ ---
    fprintf('----------------------------------------\n');
    fprintf('DAC TRUNG HINH HOC & TRUC XOAY:\n');
    fprintf('  - Ixy (Tich quan tinh):   %.2e mm4 ', Ixy);
    if abs(Ixy) > 100
        fprintf('(Tiet dien BAT DOI XUNG)\n');
    else
        fprintf('(Tiet dien DOI XUNG)\n');
    end
    fprintf('  - Goc truc chinh (Alpha): %.2f do\n', rad2deg(alpha_rad));
    fprintf('  - Goc truc trung hoa (Beta): %.2f do (So voi truc X nam ngang)\n', Beta_deg);
    fprintf('\nKET QUA UNG SUAT (GLOBAL):\n');
    fprintf('  >> MAX Sigma: %+.4f MPa\n', Real_Max);
    fprintf('  >> MIN Sigma: %+.4f MPa\n', Real_Min);
    
    if Nz ~= 0
        fprintf('  (Luu y: Do co luc doc Nz, truc trung hoa bi tinh tien, khong di qua trong tam)\n');
    end
    fprintf('----------------------------------------\n');
end
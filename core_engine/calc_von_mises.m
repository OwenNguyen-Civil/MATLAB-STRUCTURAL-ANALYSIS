function [VM_Max, SF, VM_Data] = calc_von_mises(Props, Loads, Shear_Res, Fy_Limit)
% CALC_VON_MISES (Pixel-Compatible Version)
% Tương thích với calc_shear_stress phương pháp "Máy in"
% -------------------------------------------------------------------------
% INPUT:
%   Props      : Struct đặc trưng hình học (Ix, Iy, Ixy, Area, Cx, Cy...)
%   Loads      : Tải trọng (Nz, Mx, My)
%   Shear_Res  : Kết quả từ calc_shear_stress (biến Debug_Data cũ)
%                (Phải chứa: .Y và .Tau)
%   Fy_Limit   : Giới hạn chảy của vật liệu (MPa)
%
% OUTPUT:
%   VM_Max     : Ứng suất Von Mises lớn nhất (MPa)
%   SF         : Hệ số an toàn
%   VM_Data    : Struct dữ liệu để vẽ biểu đồ
% -------------------------------------------------------------------------

    % 1. LẤY DỮ LIỆU TỪ MODULE CẮT (Mapping lại tên biến)
    % Code shear của ông trả về .Y và .Tau, ta lấy ra dùng
    Y_list   = Shear_Res.Y; 
    Tau_list = Shear_Res.Tau; 
    
    n = length(Y_list);
    
    % Khởi tạo mảng kết quả
    Sigma_Z_List  = zeros(n, 1);
    Sigma_VM_List = zeros(n, 1);
    
    % 2. CHUẨN BỊ CÔNG THỨC ỨNG SUẤT PHÁP (Navier Tổng Quát)
    Cx = Props.Cx; Cy = Props.Cy;
    Ix = Props.Ix; Iy = Props.Iy; Ixy = Props.Ixy;
    Det_I = Ix*Iy - Ixy^2;
    
    % Các hệ số của phương trình mặt phẳng ứng suất: Sigma = Const + K_y*y + K_x*x
    Const = Loads.Nz / Props.Area;
    K_y = -(Loads.Mx * Iy - Loads.My * Ixy) / Det_I; 
    K_x =  (Loads.My * Ix - Loads.Mx * Ixy) / Det_I;
    
    % 3. QUÉT DỌC TRỤC Y ĐỂ TÍNH VON MISES
    for i = 1:n
        % Lấy tọa độ Y tại lát cắt (so với trọng tâm)
        y_world = Y_list(i);
        y_local = y_world - Cy; 
        
        % Lấy Ứng suất tiếp Tau tại lát cắt này
        tau_val = Tau_list(i);
        
        % --- CHIẾN THUẬT TÍNH SIGMA ---
        % Tau Zhuravskii là giá trị trung bình trên bề rộng b(y).
        % Để kiểm tra tương tác nguy hiểm nhất, ta thường tính Sigma tại 
        % TRỤC TRỌNG TÂM CỦA LÁT CẮT (tức là tại x_local = 0).
        % Lý do: Tại mép ngoài biên, Tau thường = 0, chỉ có Sigma max.
        % Tại giữa tiết diện, Tau max, Sigma có thể nhỏ nhưng sự kết hợp mới nguy hiểm.
        
        x_local = 0; % Tính tại trục đi qua trọng tâm
        
        % Tính Sigma pháp
        sigma_val = Const + K_y * y_local + K_x * x_local;
        
        % TÍNH VON MISES: S_vm = sqrt(Sigma^2 + 3*Tau^2)
        vm_val = sqrt(sigma_val^2 + 3 * tau_val^2);
        
        % Lưu kết quả
        Sigma_Z_List(i)  = sigma_val;
        Sigma_VM_List(i) = vm_val;
    end
    
    % 4. TỔNG HỢP KẾT QUẢ
    [VM_Max, idx_max] = max(Sigma_VM_List);
    
    % Tính hệ số an toàn (Safety Factor)
    if VM_Max > 1e-6
        SF = Fy_Limit / VM_Max;
    else
        SF = 999; % Vô cùng an toàn
    end
    
    % 5. ĐÓNG GÓI DỮ LIỆU VẼ
    VM_Data.Y_coords = Y_list;
    VM_Data.Sigma    = Sigma_Z_List;
    VM_Data.Tau      = Tau_list;
    VM_Data.VonMises = Sigma_VM_List;
    VM_Data.Y_Max    = Y_list(idx_max);
    
end
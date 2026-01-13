%% MAIN SCRIPT: REPORT OPTIMIZED 
clc; clear; close all;
warning('off', 'all'); % Tắt toàn bộ cảnh báo đỏ cho đỡ rác mắt

% --- 1. INPUT DATA (GIỮ NGUYÊN) ---
Poly = [
    % --- ĐẦU & SỪNG (HEAD & HORNS) ---
    0,   120;   % Mũi nhọn
    5,   125;   % Gờ mũi trên
    15,  122;   % Hõm mũi
    20,  130;   % Gờ mắt trước
    25,  128;   % Hốc mắt
    30,  135;   % Gờ mắt sau
    35,  145;   % Sừng nhỏ 1
    42,  142;   % Chân sừng 1
    50,  158;   % Sừng chính lớn (Đỉnh đầu)
    58,  148;   % Chân sừng chính
    65,  152;   % Sừng sau gáy
    70,  145;   % Gáy

    % --- GAI LƯNG & CỔ (NECK & SPINE SPIKES) ---
    % Tạo hình răng cưa dọc sống lưng
    80,  155; 85,  150; % Gai 1
    95,  160; 100, 155; % Gai 2
    110, 165; 115, 160; % Gai 3 (Gần vai)

    % --- CÁNH PHỨC TẠP (COMPLEX WING STRUCTURE) ---
    125, 180;   % Khớp cánh trước
    135, 220;   % Xương cánh tay chính
    145, 250;   % Đỉnh ngón cánh 1 (Cao nhất)
    155, 230;   % Màng cánh võng 1
    170, 260;   % Đỉnh ngón cánh 2
    185, 235;   % Màng cánh võng 2
    205, 255;   % Đỉnh ngón cánh 3
    215, 225;   % Màng cánh võng 3
    235, 240;   % Đỉnh ngón cánh 4 (Thấp hơn)
    240, 200;   % Nách cánh sau

    % --- LƯNG DƯỚI & GAI ĐUÔI (LOWER BACK & TAIL SPIKES) ---
    260, 190; 265, 185; % Gai lưng dưới 1
    280, 185; 285, 180; % Gai lưng dưới 2
    310, 175; 315, 170; % Gai gốc đuôi
    340, 160; 345, 155; % Gai đuôi giữa
    370, 140; 375, 135; % Gai đuôi gần chóp

    % --- CHÓP ĐUÔI (TAIL TIP FIN) ---
    400, 120;   % Gốc vây đuôi trên
    420, 130;   % Đỉnh vây đuôi trên
    430, 110;   % Khe vây đuôi
    440, 100;   % Đỉnh vây đuôi chính (Điểm xa nhất X)
    425, 80;    % Đỉnh vây đuôi dưới
    400, 90;    % Gốc vây đuôi dưới

    % --- BỤNG ĐUÔI & CHÂN SAU (UNDERSIDE & HIND LEG) ---
    360, 70;    % Bụng đuôi
    320, 60;    % Hông sau
    300, 45;    % Đùi sau
    305, 25;    % Đầu gối sau
    295, 10;    % Cổ chân sau
    % Chi tiết móng chân sau chạm đất (Y=0)
    305, 0;     % Móng 1
    300, 5;     % Khe ngón
    295, 0;     % Móng 2
    290, 5;     % Khe ngón
    285, 0;     % Móng 3
    275, 20;    % Gót chân sau

    % --- BỤNG VẢY (SCALY BELLY) ---
    % Tạo hình gợn sóng cho bụng
    250, 30; 230, 25; 
    210, 32; 190, 28;
    170, 35; 150, 30;
    
    % --- CHÂN TRƯỚC (FRONT LEG) ---
    130, 40;    % Hông trước
    120, 25;    % Đùi trước
    125, 15;    % Khuỷu chân trước
    % Chi tiết móng chân trước chạm đất (Y=0)
    135, 0;     % Móng 1
    130, 5;     % Khe ngón
    125, 0;     % Móng 2
    120, 5;     % Khe ngón
    115, 0;     % Móng 3
    105, 25;    % Cổ chân trước

    % --- NGỰC, CỔ & HÀM (CHEST, THROAT & JAW) ---
    90,  50;    % Ngực
    70,  65;    % Ức
    60,  60;    % Hõm cổ dưới
    50,  85;    % Cổ họng
    40,  95;    % Góc hàm dưới
    30,  105;   % Cằm
    25,  100;   % Gai cằm nhỏ
    15,  110;   % Răng hàm dưới nhô lên
    5,   115;   % Môi dưới
    0,   120    % Về lại mũi (Khép kín)
];
CheckPoints = []; 
Ld.Nz = -2000 * 1e3;   
Ld.Mx =  450  * 1e6;   
Ld.My =  120  * 1e6;   
Vy_Input = 150 * 1e3;

% --- 2. TÍNH TOÁN (SILENT MODE) ---
% Tính Ứng suất pháp (Chặn output thừa trong hàm nếu cần)
[Sig_Res, Props] = calc_normal_stress(Poly, CheckPoints, Ld);

% Tính Ứng suất tiếp (Đo giờ)
tic;
[Tau_max, Shear_Data] = calc_shear_stress(Poly, Vy_Input, Props);
t_shear = toc;

%Kiểm tra điều kiện bền (Ứng suất tương đương)
Material_Fy = 550; 
[VM_Max, SF, VM_Res] = calc_von_mises(Props, Ld, Shear_Data, Material_Fy);

% --- 3. IN BÁO CÁO (DASHBOARD STYLE) ---
clc; % Xóa màn hình lần cuối để hiện bảng đẹp
fprintf('\n');
fprintf('==============================================================\n');
fprintf('          KẾT QUẢ TÍNH TOÁN TIẾT DIỆN (STRUCTURAL REPORT)     \n');
fprintf('==============================================================\n');

% --- BLOCK 1: HÌNH HỌC ---
fprintf(' [1] DAC TRUNG HINH HOC (GEOMETRY)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Dien tich (Area)      : %12.2f mm2\n', abs(Props.Area));
fprintf('  > Trong tam (Cx, Cy)    : (%9.2f , %9.2f) mm\n', Props.Cx, Props.Cy);
fprintf('  > Momen Quan tinh Ix    : %12.2e mm4\n', Props.Ix);
fprintf('  > Momen Quan tinh Iy    : %12.2e mm4\n', Props.Iy);
fprintf('  > Tich Quan tinh Ixy    : %12.2e mm4', Props.Ixy);
if abs(Props.Ixy) > 1000, fprintf(' -> [BAT DOI XUNG]\n'); else, fprintf(' -> [DOI XUNG]\n'); end
fprintf('  > Goc xoay truc chinh   : %12.2f do\n', Props.Alpha_Principal_deg);
fprintf('\n');

% --- BLOCK 2: TẢI TRỌNG & ỨNG SUẤT PHÁP (ĐÃ SỬA) ---
fprintf(' [2] UNG SUAT PHAP (NORMAL STRESS - SIGMA)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Tai trong dau vao     : Nz = %.0f kN | Mx = %.0f kNm | My = %.0f kNm\n', ...
        Ld.Nz/1000, Ld.Mx/1e6, Ld.My/1e6);
fprintf('  > Goc truc trung hoa NA : %.2f do (So voi truc X)\n', Props.Beta_NeutralAxis_deg);

% --- TÍNH TOÁN LẠI MAX/MIN TOÀN BỘ TIẾT DIỆN ---
% (Dùng công thức Navier tổng quát quét qua tất cả đỉnh Poly)
Cx = Props.Cx; Cy = Props.Cy;
Ix = Props.Ix; Iy = Props.Iy; Ixy = Props.Ixy;
Det_I = Ix*Iy - Ixy^2;

K_y = -(Ld.Mx * Iy - Ld.My * Ixy) / Det_I; 
K_x =  (Ld.My * Ix - Ld.Mx * Ixy) / Det_I;
Const = Ld.Nz / Props.Area;

% Quét tất cả các đỉnh của Poly
X_v = Poly(:,1) - Cx;
Y_v = Poly(:,2) - Cy;
All_Sigma = Const + K_y .* Y_v + K_x .* X_v;

S_max_Global = max(All_Sigma);
S_min_Global = min(All_Sigma);

fprintf('  > MAX SIGMA (Keo)       : %+.2f MPa  <--- [GLOBAL MAX]\n', S_max_Global);
fprintf('  > MIN SIGMA (Nen)       : %+.2f MPa  <--- [GLOBAL MIN]\n', S_min_Global);

% Nếu có CheckPoints thì in thêm
if ~isempty(Sig_Res)
    fprintf('\n  > Chi tiet tai CheckPoints:\n');
    fprintf('    | %-10s | %-15s |\n', 'Diem #', 'Ung suat (MPa)');
    fprintf('    |------------|-----------------|\n');
    for i = 1:length(Sig_Res)
        fprintf('    | %-10d | %15.4f |\n', i, Sig_Res(i));
    end
end
fprintf('\n');

% --- BLOCK 3: ỨNG SUẤT TIẾP ---
fprintf(' [3] UNG SUAT TIEP (SHEAR STRESS - TAU)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Luc cat dau vao (Vy)  : %12.2f kN\n', Vy_Input/1000);
fprintf('  > Thoi gian xu ly       : %12.4f giay\n', t_shear);
fprintf('  > TAU MAX (Zhuravskii)  : %12.2f MPa  <--- [KET QUA QUAN TRONG]\n', Tau_max);
fprintf('==============================================================\n');
fprintf('\n');

warning('off', 'MATLAB:polyshape:repairedPoly')
% --- BLOCK 4: KIỂM TRA BỀN (DESIGN CHECK) ---
fprintf(' [4] KIEM TRA BEN TONG HOP (VON MISES)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Gioi han chay Vat lieu: %.0f MPa\n', Material_Fy);
fprintf('  > MAX VON MISES Stress  : %.2f MPa\n', VM_Max);
fprintf('  > HE SO AN TOAN (SF)    : %.2f ', SF);

if SF >= 1.1
    fprintf('--> [OK - AN TOAN]\n');
elseif SF >= 1.0
    fprintf('--> [WARNING - SAT GIOI HAN]\n');
else
    fprintf('--> [DANGER - KET CAU SAP DO!!!]\n');
end
fprintf('==============================================================\n');
fprintf('\n');
%% 4. VISUALIZATION (VẼ ĐỒ THỊ)
figure(1); clf; 
set(gcf, 'Color', 'k', 'Name', 'Stress Analysis Result'); % Giao diện tối
ax = gca; hold on; axis equal; grid on;

% Vẽ Heatmap (Màu biến thiên theo ứng suất)
patch('Vertices', Poly, 'Faces', 1:size(Poly,1), ...
      'FaceVertexCData', All_Sigma, ...      % Gán dữ liệu màu là Ứng suất
      'FaceColor', 'interp', ...             % Nội suy màu giữa các điểm
      'EdgeColor', 'w', 'LineWidth', 1.5, ...
      'DisplayName', 'Phan bo Ung suat');

colormap(jet); % Dùng dải màu Jet (Xanh -> Đỏ)
c = colorbar;  % Hiện thanh chú thích màu
c.Label.String = 'Ung suat phap (MPa)';
c.Label.Color = 'w'; c.Color = 'w';
% 4.2. Vẽ Trọng tâm (Centroid)
plot(Props.Cx, Props.Cy, 'p', 'MarkerSize', 14, ...
     'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', 'DisplayName', 'Trong tam (C)');

% 4.3. Tính & Vẽ Trục trung hòa (Neutral Axis - NA)
% --- Bắt đầu khối xử lý toán học đường thẳng ---
Cx = Props.Cx; Cy = Props.Cy;
Det_I = Props.Ix * Props.Iy - Props.Ixy^2;
Const = Ld.Nz / Props.Area;
K_y = -(Ld.Mx * Props.Iy - Ld.My * Props.Ixy) / Det_I; 
K_x =  (Ld.My * Props.Ix - Ld.Mx * Props.Ixy) / Det_I;

% Xác định vùng vẽ (Viewport)
x_bnd = [min(Poly(:,1)), max(Poly(:,1))]; 
w = diff(x_bnd);
X_NA = [x_bnd(1) - w*0.5, x_bnd(2) + w*0.5]; % Kéo dài ra ngoài biên 50%

if abs(K_y) < 1e-10  % Trường hợp trục NA thẳng đứng
    if abs(K_x) > 1e-10
        val_X = Cx - Const/K_x;
        plot([val_X val_X], ylim, 'r-.', 'LineWidth', 2, 'DisplayName', 'Truc trung hoa (NA)');
    end
else                 % Trường hợp tổng quát (Ngang hoặc Xiên)
    Y_NA = Cy - (Const + K_x .* (X_NA - Cx)) ./ K_y;
    plot(X_NA, Y_NA, 'r-.', 'LineWidth', 2, 'DisplayName', 'Truc trung hoa (NA)');
end
% --- Kết thúc khối xử lý ---

% 4.4. Cấu hình hiển thị (Zoom Fix & Style)
margin = 0.3; % Lề 30%
dx = w * margin; dy = (max(Poly(:,2)) - min(Poly(:,2))) * margin;
xlim([min(Poly(:,1)) - dx, max(Poly(:,1)) + dx]);
ylim([min(Poly(:,2)) - dy, max(Poly(:,2)) + dy]);

% Trang trí trục và tiêu đề
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 10);
xlabel('X (mm)', 'FontWeight', 'bold'); ylabel('Y (mm)', 'FontWeight', 'bold');
title('BIỂU ĐỒ PHÂN BỐ ỨNG SUẤT & TRỤC TRUNG HÒA', 'Color', 'w', 'FontSize', 12);
legend('show', 'TextColor', 'w', 'Color', 'none', 'Location', 'best');
axis equal; % <--- BẮT BUỘC: Giữ đúng tỷ lệ bản đồ
grid on;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % Nền tối cho ngầu
title('BIỂU ĐỒ ỨNG SUẤT (HIGH DETAIL)', 'Color', 'w');
hold off;

%% MAIN SCRIPT: REPORT OPTIMIZED 
clc; clear; close all;
warning('off', 'all'); % Turn off all red warnings to reduce clutter

% --- 1. INPUT DATA (KEPT ORIGINAL) ---
Poly = [
    % --- HEAD & HORNS ---
    0,    120;    % Snout tip
    5,    125;    % Upper snout ridge
    15,   122;    % Snout dip
    20,   130;    % Front brow ridge
    25,   128;    % Eye socket
    30,   135;    % Rear brow ridge
    35,   145;    % Small horn 1
    42,   142;    % Horn base 1
    50,   158;    % Main large horn (Crown)
    58,   148;    % Main horn base
    65,   152;    % Nape horn
    70,   145;    % Nape

    % --- NECK & SPINE SPIKES ---
    % Create sawtooth pattern along the spine
    80,   155; 85,   150; % Spike 1
    95,   160; 100, 155; % Spike 2
    110, 165; 115, 160; % Spike 3 (Near shoulder)

    % --- COMPLEX WING STRUCTURE ---
    125, 180;    % Front wing joint
    135, 220;    % Main humerus bone
    145, 250;    % Wing digit tip 1 (Highest)
    155, 230;    % Wing membrane webbing 1
    170, 260;    % Wing digit tip 2
    185, 235;    % Wing membrane webbing 2
    205, 255;    % Wing digit tip 3
    215, 225;    % Wing membrane webbing 3
    235, 240;    % Wing digit tip 4 (Lower)
    240, 200;    % Rear underwing axilla

    % --- LOWER BACK & TAIL SPIKES ---
    260, 190; 265, 185; % Lower back spike 1
    280, 185; 285, 180; % Lower back spike 2
    310, 175; 315, 170; % Tail base spike
    340, 160; 345, 155; % Mid-tail spike
    370, 140; 375, 135; % Near-tip tail spike

    % --- TAIL TIP FIN ---
    400, 120;    % Upper tail fin base
    420, 130;    % Upper tail fin tip
    430, 110;    % Tail fin notch
    440, 100;    % Main tail fin tip (Furthest X point)
    425, 80;     % Lower tail fin tip
    400, 90;     % Lower tail fin base

    % --- UNDERSIDE & HIND LEG ---
    360, 70;     % Tail underbelly
    320, 60;     % Rear hip
    300, 45;     % Rear thigh
    305, 25;     % Rear knee
    295, 10;     % Rear ankle
    % Hind claw details touching ground (Y=0)
    305, 0;      % Claw 1
    300, 5;      % Toe gap
    295, 0;      % Claw 2
    290, 5;      % Toe gap
    285, 0;      % Claw 3
    275, 20;     % Rear heel

    % --- SCALY BELLY ---
    % Create ripple effect for belly
    250, 30; 230, 25; 
    210, 32; 190, 28;
    170, 35; 150, 30;
     
    % --- FRONT LEG ---
    130, 40;     % Front hip
    120, 25;     % Front thigh
    125, 15;     % Front elbow
    % Front claw details touching ground (Y=0)
    135, 0;      % Claw 1
    130, 5;      % Toe gap
    125, 0;      % Claw 2
    120, 5;      % Toe gap
    115, 0;      % Claw 3
    105, 25;     % Front ankle

    % --- CHEST, THROAT & JAW ---
    90,  50;     % Chest
    70,  65;     % Sternum
    60,  60;     % Lower neck dip
    50,  85;     % Throat
    40,  95;     % Jaw angle
    30,  105;    % Chin
    25,  100;    % Small chin spike
    15,  110;    % Lower tooth protrusion
    5,   115;    % Lower lip
    0,   120     % Return to nose (Closed loop)
];
CheckPoints = []; 
Ld.Nz = -2000 * 1e3;    
Ld.Mx =  450  * 1e6;    
Ld.My =  120  * 1e6;    
Vy_Input = 150 * 1e3;

% --- 2. CALCULATION (SILENT MODE) ---
% Calculate Normal Stress (Block excess output in function if needed)
[Sig_Res, Props] = calc_normal_stress(Poly, CheckPoints, Ld);

% Calculate Shear Stress (Measure time)
tic;
[Tau_max, Shear_Data] = calc_shear_stress(Poly, Vy_Input, Props);
t_shear = toc;

% Check Structural Integrity (Equivalent Stress)
Material_Fy = 550; 
[VM_Max, SF, VM_Res] = calc_von_mises(Props, Ld, Shear_Data, Material_Fy);

% --- 3. PRINT REPORT (DASHBOARD STYLE) ---
clc; % Clear screen one last time for a clean table
fprintf('\n');
fprintf('==============================================================\n');
fprintf('          STRUCTURAL REPORT (SECTION CALCULATION)             \n');
fprintf('==============================================================\n');

% --- BLOCK 1: GEOMETRY ---
fprintf(' [1] GEOMETRIC PROPERTIES\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Area                  : %12.2f mm2\n', abs(Props.Area));
fprintf('  > Centroid (Cx, Cy)     : (%9.2f , %9.2f) mm\n', Props.Cx, Props.Cy);
fprintf('  > Moment of Inertia Ix  : %12.2e mm4\n', Props.Ix);
fprintf('  > Moment of Inertia Iy  : %12.2e mm4\n', Props.Iy);
fprintf('  > Product of Inertia Ixy: %12.2e mm4', Props.Ixy);
if abs(Props.Ixy) > 1000, fprintf(' -> [ASYMMETRIC]\n'); else, fprintf(' -> [SYMMETRIC]\n'); end
fprintf('  > Principal Axis Angle  : %12.2f deg\n', Props.Alpha_Principal_deg);
fprintf('\n');

% --- BLOCK 2: NORMAL STRESS ---
fprintf(' [2] NORMAL STRESS (SIGMA)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Input Load            : Nz = %.0f kN | Mx = %.0f kNm | My = %.0f kNm\n', ...
        Ld.Nz/1000, Ld.Mx/1e6, Ld.My/1e6);
fprintf('  > Neutral Axis Angle NA : %.2f deg (vs X-axis)\n', Props.Beta_NeutralAxis_deg);

% --- RECALCULATE GLOBAL MAX/MIN FOR ENTIRE SECTION ---
% (Use generalized Navier formula sweeping all Poly vertices)
Cx = Props.Cx; Cy = Props.Cy;
Ix = Props.Ix; Iy = Props.Iy; Ixy = Props.Ixy;
Det_I = Ix*Iy - Ixy^2;

K_y = -(Ld.Mx * Iy - Ld.My * Ixy) / Det_I; 
K_x =  (Ld.My * Ix - Ld.Mx * Ixy) / Det_I;
Const = Ld.Nz / Props.Area;

% Sweep all vertices of Poly
X_v = Poly(:,1) - Cx;
Y_v = Poly(:,2) - Cy;
All_Sigma = Const + K_y .* Y_v + K_x .* X_v;

S_max_Global = max(All_Sigma);
S_min_Global = min(All_Sigma);

fprintf('  > MAX SIGMA (Tension)   : %+.2f MPa  <--- [GLOBAL MAX]\n', S_max_Global);
fprintf('  > MIN SIGMA (Compress)  : %+.2f MPa  <--- [GLOBAL MIN]\n', S_min_Global);

% If CheckPoints exist, print details
if ~isempty(Sig_Res)
    fprintf('\n  > Details at CheckPoints:\n');
    fprintf('    | %-10s | %-15s |\n', 'Point #', 'Stress (MPa)');
    fprintf('    |------------|-----------------|\n');
    for i = 1:length(Sig_Res)
        fprintf('    | %-10d | %15.4f |\n', i, Sig_Res(i));
    end
end
fprintf('\n');

% --- BLOCK 3: SHEAR STRESS ---
fprintf(' [3] SHEAR STRESS (TAU)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Input Shear Force (Vy): %12.2f kN\n', Vy_Input/1000);
fprintf('  > Processing Time       : %12.4f sec\n', t_shear);
fprintf('  > TAU MAX (Zhuravskii)  : %12.2f MPa  <--- [CRITICAL RESULT]\n', Tau_max);
fprintf('==============================================================\n');
fprintf('\n');

warning('off', 'MATLAB:polyshape:repairedPoly')
% --- BLOCK 4: DESIGN CHECK (VON MISES) ---
fprintf(' [4] COMBINED STRENGTH CHECK (VON MISES)\n');
fprintf(' -------------------------------------------------------------\n');
fprintf('  > Material Yield Limit  : %.0f MPa\n', Material_Fy);
fprintf('  > MAX VON MISES Stress  : %.2f MPa\n', VM_Max);
fprintf('  > SAFETY FACTOR (SF)    : %.2f ', SF);

if SF >= 1.1
    fprintf('--> [OK - SAFE]\n');
elseif SF >= 1.0
    fprintf('--> [WARNING - NEAR LIMIT]\n');
else
    fprintf('--> [DANGER - STRUCTURAL COLLAPSE!!!]\n');
end
fprintf('==============================================================\n');
fprintf('\n');

%% 4. VISUALIZATION (PLOTTING)
figure(1); clf; 
set(gcf, 'Color', 'k', 'Name', 'Stress Analysis Result'); % Dark interface
ax = gca; hold on; axis equal; grid on;

% Draw Heatmap (Color varies by stress)
patch('Vertices', Poly, 'Faces', 1:size(Poly,1), ...
      'FaceVertexCData', All_Sigma, ...      % Assign color data to Stress
      'FaceColor', 'interp', ...             % Interpolate color between points
      'EdgeColor', 'w', 'LineWidth', 1.5, ...
      'DisplayName', 'Stress Distribution');

colormap(jet); % Use Jet colormap (Blue -> Red)
c = colorbar;  % Show color legend
c.Label.String = 'Normal Stress (MPa)';
c.Label.Color = 'w'; c.Color = 'w';
% 4.2. Draw Centroid
plot(Props.Cx, Props.Cy, 'p', 'MarkerSize', 14, ...
     'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', 'DisplayName', 'Centroid (C)');

% 4.3. Calculate & Draw Neutral Axis (NA)
% --- Start math block ---
Cx = Props.Cx; Cy = Props.Cy;
Det_I = Props.Ix * Props.Iy - Props.Ixy^2;
Const = Ld.Nz / Props.Area;
K_y = -(Ld.Mx * Props.Iy - Ld.My * Props.Ixy) / Det_I; 
K_x =  (Ld.My * Props.Ix - Ld.Mx * Props.Ixy) / Det_I;

% Define Viewport
x_bnd = [min(Poly(:,1)), max(Poly(:,1))]; 
w = diff(x_bnd);
X_NA = [x_bnd(1) - w*0.5, x_bnd(2) + w*0.5]; % Extend 50% beyond border

if abs(K_y) < 1e-10  % Vertical NA case
    if abs(K_x) > 1e-10
        val_X = Cx - Const/K_x;
        plot([val_X val_X], ylim, 'r-.', 'LineWidth', 2, 'DisplayName', 'Neutral Axis (NA)');
    end
else                 % General case (Horizontal or Slanted)
    Y_NA = Cy - (Const + K_x .* (X_NA - Cx)) ./ K_y;
    plot(X_NA, Y_NA, 'r-.', 'LineWidth', 2, 'DisplayName', 'Neutral Axis (NA)');
end
% --- End math block ---

% 4.4. Display Configuration (Zoom Fix & Style)
margin = 0.3; % 30% Margin
dx = w * margin; dy = (max(Poly(:,2)) - min(Poly(:,2))) * margin;
xlim([min(Poly(:,1)) - dx, max(Poly(:,1)) + dx]);
ylim([min(Poly(:,2)) - dy, max(Poly(:,2)) + dy]);

% Decorate axes and title
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 10);
xlabel('X (mm)', 'FontWeight', 'bold'); ylabel('Y (mm)', 'FontWeight', 'bold');
title('STRESS DISTRIBUTION DIAGRAM & NEUTRAL AXIS', 'Color', 'w', 'FontSize', 12);
legend('show', 'TextColor', 'w', 'Color', 'none', 'Location', 'best');
axis equal; % <--- MANDATORY: Keep aspect ratio correct
grid on;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % Dark background for style
title('STRESS DIAGRAM (HIGH DETAIL)', 'Color', 'w');
hold off;

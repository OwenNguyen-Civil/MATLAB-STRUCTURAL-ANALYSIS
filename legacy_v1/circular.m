clc; clear; close all;
fprintf('=======================================================\n');
fprintf('     CIRCULAR/PIPE STRESS CALCULATION PROGRAM (V3.4)   \n');
fprintf('     (Auto-generate high-res circular coordinates)     \n');
fprintf('=======================================================\n');

% =========================================================
% 1. INPUT LOADS
% =========================================================
fprintf('\n--- STEP 1: INPUT LOADS ---\n');
Nz = input('1. Input Normal Force Nz (N) [Positive=Tension, Negative=Compression]: ');
Mx_Nm = input('2. Input Moment Mx (N.m): ');
My_Nm = input('3. Input Moment My (N.m): ');
Mx = Mx_Nm * 1000; My = My_Nm * 1000;

% =========================================================
% 2. INPUT CIRCULAR DIMENSIONS (NO COORDINATES NEEDED)
% =========================================================
fprintf('\n--- STEP 2: INPUT CIRCULAR DIMENSIONS ---\n');
D_out = input('1. Input OUTER DIAMETER (mm): ');
D_in  = input('2. Input INNER DIAMETER (mm) [Enter 0 for Solid Shaft]: ');

% Auto-generate coordinates (100 points for smoothness)
N_points = 100; 
theta = linspace(0, 2*pi, N_points);

% Create outer circle coordinates
R_out = D_out / 2;
x_out = R_out * cos(theta);
y_out = R_out * sin(theta);

% Create inner circle coordinates (if hollow)
has_hole = false;
if D_in > 0
    has_hole = true;
    R_in = D_in / 2;
    x_in = R_in * cos(theta);
    y_in = R_in * sin(theta);
end

% =========================================================
% 3. INPUT POINTS OF INTEREST
% =========================================================
fprintf('\n--- STEP 3: SELECT POINTS OF INTEREST ---\n');
fprintf('Hint: Enter [0 R_out; 0 -R_out] to check top and bottom extreme points.\n');
pts_matrix = input('>> Enter points matrix [x y]: ');
try
    x_pts = pts_matrix(:, 1)'; y_pts = pts_matrix(:, 2)';
catch
    error('Points format error!');
end

% =========================================================
% 4. GEOMETRY PROCESSING (AUTO MESH)
% =========================================================
pgon_main = polyshape(x_out, y_out);

if has_hole
    pgon_hole = polyshape(x_in, y_in);
    pgon = subtract(pgon_main, pgon_hole);
    fprintf('\n[INFO] Created HOLLOW PIPE section.\n');
else
    pgon = pgon_main;
    fprintf('\n[INFO] Created SOLID SHAFT section.\n');
end

% Calculate geometric properties
Area = area(pgon);          
[Cx, Cy] = centroid(pgon);  

% Calculate I (Exact formula for circle for comparison)
% However, we still use the polygon algorithm to match Heatmap data
% This inertia calculation uses the polygon algorithm (like in hollow.m)
% -----------------------------------------------------------
[Cx_out, Cy_out] = centroid(pgon_main); Area_out = area(pgon_main);
Ixx_out = 0; Iyy_out = 0;
x_loop = [x_out, x_out(1)]; y_loop = [y_out, y_out(1)];
for i = 1:length(x_out)
    Xi = x_loop(i); Xi1 = x_loop(i+1); Yi = y_loop(i); Yi1 = y_loop(i+1);
    D = Xi*Yi1 - Xi1*Yi;
    Ixx_out = Ixx_out + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
    Iyy_out = Iyy_out + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
end
dx_out = Cx_out - Cx; dy_out = Cy_out - Cy;
Ix_total = (Ixx_out - Area_out*Cy_out^2) + Area_out*dy_out^2;
Iy_total = (Iyy_out - Area_out*Cx_out^2) + Area_out*dx_out^2;

if has_hole
    [Cx_in, Cy_in] = centroid(pgon_hole); Area_in = area(pgon_hole);
    Ixx_in = 0; Iyy_in = 0;
    x_loop_in = [x_in, x_in(1)]; y_loop_in = [y_in, y_in(1)];
    for i = 1:length(x_in)
        Xi = x_loop_in(i); Xi1 = x_loop_in(i+1); Yi = y_loop_in(i); Yi1 = y_loop_in(i+1);
        D = Xi*Yi1 - Xi1*Yi;
        Ixx_in = Ixx_in + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
        Iyy_in = Iyy_in + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
    end
    dx_in = Cx_in - Cx; dy_in = Cy_in - Cy;
    Ix_hole_shifted = (Ixx_in - Area_in*Cy_in^2) + Area_in*dy_in^2;
    Iy_hole_shifted = (Iyy_in - Area_in*Cx_in^2) + Area_in*dx_in^2;
    
    Ix = Ix_total - Ix_hole_shifted; Iy = Iy_total - Iy_hole_shifted;
else
    Ix = Ixx_out - Area_out * Cy_out^2; Iy = Iyy_out - Area_out * Cx_out^2;
end
% -----------------------------------------------------------

% =========================================================
% 5. RESULT REPORT
% =========================================================
fprintf('\n================ REPORT ================\n');
fprintf('1. Geometric Properties:\n');
fprintf('   - Area A:            %.2f mm^2\n', Area);
fprintf('   - Ix (Horiz. Axis):  %.3e mm^4\n', Ix);
fprintf('   - Iy (Vert. Axis):   %.3e mm^4\n', Iy);

% Quick check with theoretical formula (pi*D^4/64)
I_theory = (pi*(D_out^4 - D_in^4))/64;
fprintf('   (Theory Check:       %.3e mm^4) -> Error: %.4f%%\n', I_theory, abs(Ix-I_theory)/I_theory*100);
fprintf('----------------------------------------\n');

Points_X_new = x_pts - Cx; Points_Y_new = y_pts - Cy;
Stress_Specific = zeros(1, length(x_pts)); 
for k = 1:length(x_pts)
    val = (Nz / Area) + (Mx * Points_Y_new(k) / Ix) + (My * Points_X_new(k) / Iy);
    Stress_Specific(k) = val;
    fprintf('Point %d (X=%.1f, Y=%.1f) | Sigma = %8.3f MPa\n', k, x_pts(k), y_pts(k), val);
end

% =========================================================
% 6. DRAW CIRCULAR HEATMAP (USING PCOLOR)
% =========================================================
figure('Name', 'Circular Stress Map', 'Color', 'k');
pgon_shifted = translate(pgon, -Cx, -Cy);
[xlim, ylim] = boundingbox(pgon_shifted);
res = 200; 
xg = linspace(xlim(1), xlim(2), res); yg = linspace(ylim(1), ylim(2), res);
[X_grid, Y_grid] = meshgrid(xg, yg);
X_vec = X_grid(:); Y_vec = Y_grid(:);

in_vec = isinterior(pgon_shifted, X_vec, Y_vec);
in = reshape(in_vec, size(X_grid));
Stress_Grid = nan(size(X_grid)); 
Stress_Grid(in) = (Nz / Area) + (Mx .* Y_grid(in) ./ Ix) + (My .* X_grid(in) ./ Iy);

h = pcolor(X_grid, Y_grid, Stress_Grid);
set(h, 'EdgeColor', 'none'); shading interp; hold on; 

plot(pgon_shifted, 'FaceColor', 'none', 'EdgeColor', 'w', 'LineWidth', 2);
xline(0, 'w-.'); yline(0, 'w-.');
scatter(Points_X_new, Points_Y_new, 80, 'w', 'filled', 'MarkerEdgeColor', 'k');

colormap(jet); c = colorbar; c.Color = 'w'; c.Label.Color = 'w';
c.Label.String = 'Stress (MPa)';
title('CIRCULAR/PIPE STRESS DIAGRAM', 'Color', 'w'); axis equal; grid on;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');

for k = 1:length(x_pts)
    text(Points_X_new(k), Points_Y_new(k), sprintf(' %.2f', Stress_Specific(k)), ...
        'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'w');
end
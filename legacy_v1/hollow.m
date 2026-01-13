clc; clear; close all;
fprintf('=======================================================\n');
fprintf('      GENERAL STRESS CALCULATION PROGRAM (V3.3)        \n');
fprintf('     (New Feature: Supports HOLLOW/PERFORATED Sections) \n');
fprintf('=======================================================\n');

% =========================================================
% 1. INPUT LOADS
% =========================================================
fprintf('\n--- STEP 1: INPUT LOADS ---\n');
Nz = input('1. Input Normal Force Nz (N) [Positive=Tension, Negative=Compression]: ');
Mx_Nm = input('2. Input Moment Mx (N.m) [Bending about Horizontal Axis]: ');
My_Nm = input('3. Input Moment My (N.m) [Bending about Vertical Axis]: ');

Mx = Mx_Nm * 1000; % N.mm
My = My_Nm * 1000; % N.mm

% =========================================================
% 2. INPUT SECTION SHAPE (HOLLOW/HOLE SUPPORTED)
% =========================================================
fprintf('\n--- STEP 2: INPUT SECTION SHAPE ---\n');

% A. Input Outer Boundary
fprintf('A. Input OUTER BOUNDARY coordinates [x y] (e.g., Rect 100x200):\n');
fprintf('   [0 0; 100 0; 100 200; 0 200]\n');
outer_matrix = input('>> Outer boundary matrix: ');
try
    x_out = outer_matrix(:, 1)'; 
    y_out = outer_matrix(:, 2)';
catch
    error('Outer boundary format error!');
end

% B. Input Hole Coordinates
fprintf('\nB. Input HOLE COORDINATES [x y] (Press Enter to skip if solid):\n');
fprintf('   Example hole in center: [20 20; 80 20; 80 180; 20 180]\n');
inner_matrix = input('>> Hole matrix: ');

has_hole = false;
if ~isempty(inner_matrix)
    try
        x_in = inner_matrix(:, 1)'; 
        y_in = inner_matrix(:, 2)';
        has_hole = true;
    catch
        error('Hole format error!');
    end
end

% =========================================================
% 3. INPUT POINTS OF INTEREST
% =========================================================
fprintf('\n--- STEP 3: INPUT POINTS OF INTEREST [x y] ---\n');
pts_matrix = input('>> Enter points of interest matrix: ');

try
    x_pts = pts_matrix(:, 1)';
    y_pts = pts_matrix(:, 2)';
catch
    error('Points of interest format error!');
end

% =========================================================
% 4. GEOMETRY & PROPERTIES PROCESSING (CORE ENGINE)
% =========================================================
% Create geometric objects
pgon_main = polyshape(x_out, y_out);

if has_hole
    pgon_hole = polyshape(x_in, y_in);
    pgon = subtract(pgon_main, pgon_hole); % Geometric subtraction
    fprintf('\n[INFO] Processing HOLLOW section...\n');
else
    pgon = pgon_main;
    fprintf('\n[INFO] Processing SOLID section...\n');
end

% Calculate actual Area and Centroid
Area = area(pgon);          
[Cx, Cy] = centroid(pgon);  

% --- MOMENT OF INERTIA ALGORITHM (STEINER FOR HOLLOW SHAPES) ---
% I_actual = I_outer - I_hole (all shifted to common Centroid C)

% 1. Calculate for Outer Boundary
[Cx_out, Cy_out] = centroid(pgon_main);
Area_out = area(pgon_main);
Ixx_out = 0; Iyy_out = 0;
x_loop = [x_out, x_out(1)]; y_loop = [y_out, y_out(1)];
for i = 1:length(x_out)
    Xi = x_loop(i); Xi1 = x_loop(i+1); Yi = y_loop(i); Yi1 = y_loop(i+1);
    D = Xi*Yi1 - Xi1*Yi;
    Ixx_out = Ixx_out + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
    Iyy_out = Iyy_out + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
end
% Shift to common Centroid
dx_out = Cx_out - Cx; dy_out = Cy_out - Cy;
Ix_total = (Ixx_out - Area_out*Cy_out^2) + Area_out*dy_out^2;
Iy_total = (Iyy_out - Area_out*Cx_out^2) + Area_out*dx_out^2;

% 2. Calculate for Hole (if any) and subtract
if has_hole
    [Cx_in, Cy_in] = centroid(pgon_hole);
    Area_in = area(pgon_hole);
    Ixx_in = 0; Iyy_in = 0;
    x_loop_in = [x_in, x_in(1)]; y_loop_in = [y_in, y_in(1)];
    for i = 1:length(x_in)
        Xi = x_loop_in(i); Xi1 = x_loop_in(i+1); Yi = y_loop_in(i); Yi1 = y_loop_in(i+1);
        D = Xi*Yi1 - Xi1*Yi;
        Ixx_in = Ixx_in + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
        Iyy_in = Iyy_in + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
    end
    % Shift to common Centroid
    dx_in = Cx_in - Cx; dy_in = Cy_in - Cy;
    Ix_hole_shifted = (Ixx_in - Area_in*Cy_in^2) + Area_in*dy_in^2;
    Iy_hole_shifted = (Iyy_in - Area_in*Cx_in^2) + Area_in*dx_in^2;
    
    % FINAL RESULTS
    Ix = Ix_total - Ix_hole_shifted;
    Iy = Iy_total - Iy_hole_shifted;
else
    Ix = Ixx_out - Area_out * Cy_out^2;
    Iy = Iyy_out - Area_out * Cx_out^2;
end

% =========================================================
% 5. RESULT REPORT
% =========================================================
fprintf('\n=================================================\n');
fprintf('             RESULT REPORT                       \n');
fprintf('=================================================\n');
fprintf('1. Geometric Properties:\n');
fprintf('   - Area A:            %.2f mm^2\n', Area);
fprintf('   - Centroid C:        (%.2f, %.2f) mm\n', Cx, Cy);
fprintf('   - Ix (Horiz. Axis):  %.3e mm^4\n', Ix);
fprintf('   - Iy (Vert. Axis):   %.3e mm^4\n', Iy);
fprintf('-------------------------------------------------\n');
fprintf('2. Stress at Specific Points (MPa):\n');

Points_X_new = x_pts - Cx;
Points_Y_new = y_pts - Cy;
Stress_Specific = zeros(1, length(x_pts)); 

for k = 1:length(x_pts)
    x = Points_X_new(k);
    y = Points_Y_new(k);
    val = (Nz / Area) + (Mx * y / Ix) + (My * x / Iy);
    Stress_Specific(k) = val;
    fprintf('Point %d (X=%.1f, Y=%.1f) | Sigma = %8.3f MPa\n', k, x_pts(k), y_pts(k), val);
end

% =========================================================
% 6. DRAW HEATMAP (BLACK SCREEN FIX - USING PCOLOR)
% =========================================================
figure('Name', 'Stress Contour Map (Hollow Fix)', 'Color', 'k');

% a. Process plotting data
pgon_shifted = translate(pgon, -Cx, -Cy);
[xlim, ylim] = boundingbox(pgon_shifted);

% Increase resolution for smoothness (300)
res = 300; 
xg = linspace(xlim(1), xlim(2), res);
yg = linspace(ylim(1), ylim(2), res);
[X_grid, Y_grid] = meshgrid(xg, yg);
X_vec = X_grid(:); Y_vec = Y_grid(:);

% Check interior points (isinterior handles holes correctly)
in_vec = isinterior(pgon_shifted, X_vec, Y_vec);
in = reshape(in_vec, size(X_grid));

% Calculate stress
Stress_Grid = nan(size(X_grid)); 
% Calculate only for solid parts (in), hollow parts remain NaN
Stress_Grid(in) = (Nz / Area) + (Mx .* Y_grid(in) ./ Ix) + (My .* X_grid(in) ./ Iy);

% b. Draw using PCOLOR (Better than contourf for hollow shapes)
h = pcolor(X_grid, Y_grid, Stress_Grid);
set(h, 'EdgeColor', 'none'); % Remove grid edges for smooth color
shading interp; % Interpolate colors
hold on; 

% c. Draw white section boundary (Including holes)
plot(pgon_shifted, 'FaceColor', 'none', 'EdgeColor', 'w', 'LineWidth', 2);

% d. Neutral Axis & Points of Interest
xline(0, 'w-.', 'LineWidth', 1); yline(0, 'w-.', 'LineWidth', 1);
scatter(Points_X_new, Points_Y_new, 80, 'w', 'filled', 'MarkerEdgeColor', 'k');

% e. Dark Mode Decoration
colormap(jet); 
c = colorbar; c.Label.String = 'Stress (MPa)'; 
c.Color = 'w'; c.Label.Color = 'w';

% Handle if all stress = 0 (Avoid black screen when no load)
if max(abs(Stress_Grid(:)), [], 'omitnan') == 0
    caxis([-1 1]); % Create dummy range to show green color (0)
    title('DIAGRAM (NO LOAD)', 'Color', 'w');
else
    title('STRESS DISTRIBUTION DIAGRAM (HEATMAP)', 'Color', 'w');
end

xlabel('mm', 'Color', 'w'); ylabel('mm', 'Color', 'w'); 
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3); 
axis equal; grid on;

% f. Display values
for k = 1:length(x_pts)
    text(Points_X_new(k), Points_Y_new(k), sprintf('  %.2f', Stress_Specific(k)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'w');
end
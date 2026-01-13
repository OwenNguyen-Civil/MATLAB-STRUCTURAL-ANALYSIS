clc; clear; close all;
fprintf('=======================================================\n');
fprintf('      GENERAL STRESS CALCULATION PROGRAM (V3.2)        \n');
% =========================================================
% 1. INPUT LOADS
% =========================================================
fprintf('\n--- STEP 1: INPUT LOADS ---\n');
Nz = input('1. Input Normal Force Nz (N) [Positive=Tension, Negative=Compression]: ');
Mx_Nm = input('2. Input Moment Mx (N.m) [Bending about Horizontal Axis]: ');
My_Nm = input('3. Input Moment My (N.m) [Bending about Vertical Axis]: ');

Mx = Mx_Nm * 1000; % Convert to N.mm
My = My_Nm * 1000; % Convert to N.mm

% =========================================================
% 2. INPUT VERTEX COORDINATES
% =========================================================
fprintf('\n--- STEP 2: INPUT SECTION VERTEX COORDINATES [x y] ---\n');
fprintf('Example: [0 0; 100 0; 100 50; 0 50]\n');
poly_matrix = input('>> Enter vertex coordinates matrix: ');

try
    x_poly = poly_matrix(:, 1)'; 
    y_poly = poly_matrix(:, 2)'; 
catch
    error('FORMAT ERROR: Please enter a valid 2-column matrix [x y]!');
end

% =========================================================
% 3. INPUT SPECIFIC POINTS FOR CALCULATION
% =========================================================
fprintf('\n--- STEP 3: INPUT POINTS OF INTEREST [x y] ---\n');
pts_matrix = input('>> Enter specific points matrix: ');

try
    x_pts = pts_matrix(:, 1)';
    y_pts = pts_matrix(:, 2)';
catch
    error('FORMAT ERROR: Please enter a valid 2-column matrix [x y]!');
end

% =========================================================
% 4. GEOMETRY PROCESSING (CORE ENGINE)
% =========================================================
pgon = polyshape(x_poly, y_poly);
Area = area(pgon);          
[Cx, Cy] = centroid(pgon);  

Ixx_O = 0; Iyy_O = 0; n = length(x_poly);
if x_poly(1) ~= x_poly(end) || y_poly(1) ~= y_poly(end)
    x_loop = [x_poly, x_poly(1)]; y_loop = [y_poly, y_poly(1)]; n_edge = n;
else
    x_loop = x_poly; y_loop = y_poly; n_edge = n - 1;
end

for i = 1:n_edge
    Xi = x_loop(i); Xi1 = x_loop(i+1); Yi = y_loop(i); Yi1 = y_loop(i+1);
    D = Xi*Yi1 - Xi1*Yi;
    Ixx_O = Ixx_O + (Yi^2 + Yi*Yi1 + Yi1^2)*D/12;
    Iyy_O = Iyy_O + (Xi^2 + Xi*Xi1 + Xi1^2)*D/12;
end
Ix = Ixx_O - Area * Cy^2; 
Iy = Iyy_O - Area * Cx^2; 

% =========================================================
% 5. CALCULATION & RESULT REPORTING
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
% 6. DRAW HEATMAP (DARK MODE)
% =========================================================
% Set figure background color to black ('k')
figure('Name', 'Stress Contour Map (Dark)', 'Color', 'k');

% a. Prepare plotting data
pgon_shifted = translate(pgon, -Cx, -Cy);
[xlim, ylim] = boundingbox(pgon_shifted);
res = 200; 
xg = linspace(xlim(1), xlim(2), res);
yg = linspace(ylim(1), ylim(2), res);
[X_grid, Y_grid] = meshgrid(xg, yg);
X_vec = X_grid(:); Y_vec = Y_grid(:);
in_vec = isinterior(pgon_shifted, X_vec, Y_vec);
in = reshape(in_vec, size(X_grid));
Stress_Grid = nan(size(X_grid)); 
Stress_Grid(in) = (Nz / Area) + (Mx .* Y_grid(in) ./ Ix) + (My .* X_grid(in) ./ Iy);

% b. Draw Contour
contourf(X_grid, Y_grid, Stress_Grid, 50, 'LineStyle', 'none'); 
hold on; 

% c. Draw section boundary in WHITE ('w')
plot(pgon_shifted, 'FaceColor', 'none', 'EdgeColor', 'w', 'LineWidth', 2);

% d. Draw neutral axis in WHITE
xline(0, 'w-.', 'LineWidth', 1); yline(0, 'w-.', 'LineWidth', 1);

% e. Mark specific points (White dot with black edge)
scatter(Points_X_new, Points_Y_new, 80, 'w', 'filled', 'MarkerEdgeColor', 'k');

% f. Dark Mode Decoration
colormap(jet); 
c = colorbar; 
c.Label.String = 'Stress (MPa)'; 
c.Color = 'w'; % Text on colorbar in white
c.Label.Color = 'w';

title('STRESS DISTRIBUTION DIAGRAM (HEATMAP)', 'Color', 'w');
xlabel('mm', 'Color', 'w'); ylabel('mm', 'Color', 'w'); 
axis equal; grid on;

% Set colors for axes and grid (White on Black)
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);

% g. Display values (Black background, white text)
for k = 1:length(x_pts)
    text(Points_X_new(k), Points_Y_new(k), sprintf('  %.2f', Stress_Specific(k)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'w');
end
function [Tau_max, Debug_Data] = calc_shear_stress(PolyCoords, Vy, Geom_Props)
% CALC_SHEAR_STRESS: Calculate shear stress using the "Printer" method (Pixelization)
% -----------------------------------------------------------------------
% INPUTS:
%   PolyCoords : Matrix [x y] of section vertices
%   Vy         : Vertical shear force (N)
%   Geom_Props : Struct containing geometric properties (from previous function)
%                (Requires: Cx, Cy, Ix)
%
% OUTPUTS:
%   Tau_max    : Maximum shear stress (MPa)
%   Debug_Data : Struct containing data for plotting (Y_grid, Tau_profile)
% -----------------------------------------------------------------------

    % 1. "PRINTER" CONFIGURATION (RESOLUTION)
    res_y = 500; % Resolution of 500 vertical slices
    y_poly = PolyCoords(:, 2);
    x_poly = PolyCoords(:, 1);
    
    min_y = min(y_poly); 
    max_y = max(y_poly);
    
    % Create y-slices (from bottom to top)
    y_slices = linspace(min_y, max_y, res_y)';
    dy = y_slices(2) - y_slices(1);
    
    % Create horizontal grid (to count width b)
    % Extend slightly beyond bounds to ensure the shape is fully covered
    min_x = min(x_poly); max_x = max(x_poly);
    x_grid_line = linspace(min_x, max_x, 200); 
    dx = x_grid_line(2) - x_grid_line(1);
    
    [X_grid, Y_grid] = meshgrid(x_grid_line, y_slices);
    
    % 2. SCAN SHAPE (PIXELIZATION)
    pgon = polyshape(x_poly, y_poly);
    in = isinterior(pgon, X_grid(:), Y_grid(:));
    
    % Convert 'in' vector to a matrix matching grid dimensions
    Mask = reshape(in, size(X_grid));
    
    % 3. MECHANICAL CALCULATION (PHYSICS LOOP)
    % a. Calculate cut width b(y) at each height
    % b(y) = Number of horizontal pixels * dx
    b_profile = sum(Mask, 2) * dx;
    
    % b. Calculate First Moment of Area Sx(y)
    % Sx = Integral (y - Cy) * dA from y to top edge
    % dA of a row = b_profile * dy
    
    Y_dist = y_slices - Geom_Props.Cy; % Distance to Neutral Axis
    Area_row = b_profile * dy;         % Area of each thin slice
    
    dSx = Y_dist .* Area_row;          % First moment of area of each slice
    
    % Accumulate from top down (Flip -> Cumsum -> Flip)
    Sx_profile = flipud(cumsum(flipud(dSx)));
    
    % 4. CALCULATE SHEAR STRESS (FORMULA: Tau = V*Sx / I*b)
    Tau_profile = zeros(size(Sx_profile));
    
    % Calculate only where material exists (b > 0) to avoid division by zero
    valid_idx = b_profile > 1e-6; 
    
    % Zhuravskii Formula
    Tau_profile(valid_idx) = (abs(Vy) * abs(Sx_profile(valid_idx))) ./ ...
                             (Geom_Props.Ix * b_profile(valid_idx));
                             
    % 5. RESULTS
    Tau_max = max(Tau_profile);
    
    % Pack data for plotting (if needed)
    Debug_Data.Y = y_slices;
    Debug_Data.Sx = Sx_profile;
    Debug_Data.b = b_profile;
    Debug_Data.Tau = Tau_profile;
end

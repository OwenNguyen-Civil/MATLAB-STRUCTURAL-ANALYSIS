function [VM_Max, SF, VM_Data] = calc_von_mises(Props, Loads, Shear_Res, Fy_Limit)
% CALC_VON_MISES (Pixel-Compatible Version)
% Compatible with the "Printer" method in calc_shear_stress
% -------------------------------------------------------------------------
% INPUT:
%   Props      : Geometric properties struct (Ix, Iy, Ixy, Area, Cx, Cy...)
%   Loads      : Applied Loads (Nz, Mx, My)
%   Shear_Res  : Results from calc_shear_stress (formerly Debug_Data)
%                (Must contain: .Y and .Tau)
%   Fy_Limit   : Material Yield Strength (MPa)
%
% OUTPUT:
%   VM_Max     : Maximum Von Mises Stress (MPa)
%   SF         : Safety Factor
%   VM_Data    : Struct containing data for plotting
% -------------------------------------------------------------------------

    % 1. RETRIEVE DATA FROM SHEAR MODULE (Remap variables)
    % The shear function returns .Y and .Tau, extract them for use
    Y_list   = Shear_Res.Y; 
    Tau_list = Shear_Res.Tau; 
    
    n = length(Y_list);
    
    % Initialize result arrays
    Sigma_Z_List  = zeros(n, 1);
    Sigma_VM_List = zeros(n, 1);
    
    % 2. PREPARE NORMAL STRESS FORMULA (Generalized Navier)
    Cx = Props.Cx; Cy = Props.Cy;
    Ix = Props.Ix; Iy = Props.Iy; Ixy = Props.Ixy;
    Det_I = Ix*Iy - Ixy^2;
    
    % Coefficients of the stress plane equation: Sigma = Const + K_y*y + K_x*x
    Const = Loads.Nz / Props.Area;
    K_y = -(Loads.Mx * Iy - Loads.My * Ixy) / Det_I; 
    K_x =  (Loads.My * Ix - Loads.Mx * Ixy) / Det_I;
    
    % 3. SWEEP ALONG Y-AXIS TO CALCULATE VON MISES
    for i = 1:n
        % Get Y coordinate at slice (relative to centroid)
        y_world = Y_list(i);
        y_local = y_world - Cy; 
        
        % Get Shear Stress Tau at this slice
        tau_val = Tau_list(i);
        
        % --- SIGMA CALCULATION STRATEGY ---
        % Zhuravskii Tau is the average value over the width b(y).
        % To check the most dangerous interaction, we typically calculate Sigma 
        % at the CENTROIDAL AXIS OF THE SLICE (i.e., at x_local = 0).
        % Reason: At the outer edges, Tau is usually 0 (only max Sigma exists).
        % At the center of the section, Tau is max; Sigma might be small, 
        % but the combination is what makes it critical.
        
        x_local = 0; % Calculate at the axis passing through the centroid
        
        % Calculate Normal Sigma
        sigma_val = Const + K_y * y_local + K_x * x_local;
        
        % CALCULATE VON MISES: S_vm = sqrt(Sigma^2 + 3*Tau^2)
        vm_val = sqrt(sigma_val^2 + 3 * tau_val^2);
        
        % Store results
        Sigma_Z_List(i)  = sigma_val;
        Sigma_VM_List(i) = vm_val;
    end
    
    % 4. SUMMARIZE RESULTS
    [VM_Max, idx_max] = max(Sigma_VM_List);
    
    % Calculate Safety Factor (SF)
    if VM_Max > 1e-6
        SF = Fy_Limit / VM_Max;
    else
        SF = 999; % Infinitely safe
    end
    
    % 5. PACK PLOTTING DATA
    VM_Data.Y_coords = Y_list;
    VM_Data.Sigma    = Sigma_Z_List;
    VM_Data.Tau      = Tau_list;
    VM_Data.VonMises = Sigma_VM_List;
    VM_Data.Y_Max    = Y_list(idx_max);
    
end

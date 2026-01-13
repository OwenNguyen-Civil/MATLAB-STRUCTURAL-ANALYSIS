% --- 1. SETUP ---
addpath('Core engine'); % Load the library

% Define Geometry (Coordinates of the Polygon)
Poly = [400,300; 0,300; 0,250; 80,250; 80,50; ...
        100,50; 100,250; 300,250; 300,50; ...
        320,50; 320,250; 400,250];

% Define Loads
Loads.Nz = -500 * 1e3;   % Axial Compression (N)
Loads.Mx =  30  * 1e6;   % Moment X (N.mm)
Loads.My =  20  * 1e6;   % Moment Y (N.mm)
Vy_Input =  60  * 1e3;   % Shear Force (N)
Fy_Limit =  345;         % Yield Strength (MPa)

% --- 2. SOLVE ---
% Calculate Normal Stress & Geometry
[Sig_Res, Props] = calc_normal_stress(Poly, [], Loads);

% Calculate Shear Stress (Zhuravskii)
[Tau_max, Shear_Data] = calc_shear_stress(Poly, Vy_Input, Props);

% Check Safety Factor (Von Mises)
[VM_Max, SF, VM_Res] = calc_von_mises(Props, Loads, Shear_Data, Fy_Limit);

% --- 3. RESULT ---
fprintf('Max Von Mises: %.2f MPa\n', VM_Max);
fprintf('Safety Factor: %.2f\n', SF);
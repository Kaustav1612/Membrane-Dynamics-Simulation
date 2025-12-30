% Membrane Dynamics Simulation with Physical Mode Evolution and Advanced Visualization
clear;close all; clc;

%% ================= PHYSICAL PARAMETERS =================

simulation_grid_size = 60;       % Reduced for faster testing
simulated_area =5e-12;        % 5μm x 5μm membrane (15 μm²)
membrane_area = 5e-12;            % Membrane area should match simulated area
viscosity = 140;                 % Water viscosity (Pa·s)
temperature = 310;                % Temperature (K)
boltzmann_const = 1.38e-23;       % Boltzmann constant (J/K)
bending_rigidity = 16 * boltzmann_const * temperature; % ~20 kT
time_step = 10e-6;                 % MUCH smaller time step (20 ms)
total_simulation_time = 5;     % 5 s total timefTpo
surface_tension = 60e-6;           % MUCH lower tension (0.1 pN/um)
mu = 1e-6;                       % Shear Modulus (N/m)                  
old_dir = pwd;
dir = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_absent\homogenous\01(non-sine)';
output_dir_irm = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_absent\homogenous\01(non-sine)\results_irm';            
output_dir_moransI = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_absent\homogenous\01(non-sine)\results_moransI';
output_dir_curvature_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_absent\homogenous\01(non-sine)\results_curvature_map';
pinning_strength = 30e-6;          % Moderate pinning (N/m²)
area_compressibility = 1e-6;    
refractive_indices = [1.33, 1.48, 1.35];  % [air,media,membrane,cytosol]
separation_distance = 5e-9;
base_intensity = 0.1;
wavelength = 546e-9;
V_ext = ones(simulation_grid_size)*1e-8;  % Example: external potential map
P_ext = ones(simulation_grid_size)*0.1e-6;                         % External pressure term
low_kappa = 15;
high_kappa = 20;
low_tension = 35e-6;
high_tension  = 75e-6;


%% ================= SPATIAL SETUP =================
system_length = sqrt(simulated_area);  % Length of membrane side (m)
grid_spacing = system_length / simulation_grid_size; % Correct grid spacing

x_coords = (0:simulation_grid_size-1) * grid_spacing;
y_coords = (0:simulation_grid_size-1) * grid_spacing;
[grid_X, grid_Y] = meshgrid(x_coords, y_coords);

% Convert to nanometers for visualization
nm_scale = 1e9;
X_nm = grid_X * nm_scale;
Y_nm = grid_Y * nm_scale;
L_nm = system_length* nm_scale;


fprintf('Grid spacing: %.1f nm\n', grid_spacing*1e9);
fprintf('Grid points: %d x %d\n', simulation_grid_size, simulation_grid_size);

%% ================= GENERATE HETEROGENEITY WITH PERLIN NOISE =================
% Define the grid of regions
num_region_rows = 100;
num_region_cols = 100;

% Calculate region size in pixels
region_height = floor(simulation_grid_size / num_region_rows);
region_width  = floor(simulation_grid_size / num_region_cols);

% Base values
base_kappa = bending_rigidity;
base_gamma = surface_tension;
kappa_map = base_kappa;
gamma_map = base_gamma;


%% ================= VISUALIZATION =================
figure();
% Plot kappa map
imagesc(1:simulation_grid_size, 1:simulation_grid_size, kappa_map/(boltzmann_const*temperature));
colormap jet;
axis square; colorbar; 
caxis([low_kappa high_kappa]);
title('Bending Rigidity (κ) Map in kBT');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Add grid lines to show region boundaries
% hold on;
% for i = 1:num_region_rows
%     y_pos = i * region_height + 0.5;
%     plot([0.5, simulation_grid_size+0.5], [y_pos, y_pos], 'r-', 'LineWidth', 1.5);
% end
% for j = 1:num_region_cols
%     x_pos = j * region_width + 0.5;
%     plot([x_pos, x_pos], [0.5, simulation_grid_size+0.5], 'r-', 'LineWidth', 1.5);
% end
% hold off;

% Plot gamma map
figure();
imagesc(1:simulation_grid_size, 1:simulation_grid_size, gamma_map*1e6);
colormap jet;
caxis([low_tension*1e6 high_tension*1e6]);
axis square; colorbar;
title('Surface Tension (γ) Map in pN/um');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Add grid lines
% hold on;
% for i = 1:num_region_rows
%     y_pos = i * region_height + 0.5;
%     plot([0.5, simulation_grid_size+0.5], [y_pos, y_pos], 'r-', 'LineWidth', 1.5);
% end
% 
% for j = 1:num_region_cols
%     x_pos = j * region_width + 0.5;
%     plot([x_pos, x_pos], [0.5, simulation_grid_size+0.5], 'r-', 'LineWidth', 1.5);
% end
hold off;
% close;

%% ================= MEMBRANE PINNING =================
pinning_grid_spacing = 1;  % Pin every 20th grid point
pinned_sites = false(simulation_grid_size, simulation_grid_size);
pinned_sites(1:pinning_grid_spacing:end, 1:pinning_grid_spacing:end) = true;

% Soft harmonic confinement at edges

[X_grid, Y_grid] = meshgrid(1:simulation_grid_size, 1:simulation_grid_size);

% Combine pinning and confinement
pinning_potential = zeros(simulation_grid_size, simulation_grid_size);
pinning_potential(pinned_sites) = pinning_strength;

%% ================= TIME EVOLUTION SETUP =================
num_time_steps = round(total_simulation_time / time_step);
time_points = (0:num_time_steps-1) * time_step;


% Define the amplitude of the sine wave (e.g., 5 nm)
A = 5e-9; 

% --- MODIFIED PART ---
% Desired Wavelength (lambda) is 10 times the grid size (L)
% Calculate the corresponding wave number (k = 2*pi / lambda)
k = (2 * pi) / (10 * grid_spacing); 
% --- END MODIFIED PART ---

% Generate the spatial coordinates (x and y) for the grid
% Coordinates run from 1 to L
[X, Y] = meshgrid(1:simulation_grid_size, 1:simulation_grid_size);

% Generate a single random phase offset (theta_rand) between 0 and 2*pi
theta_rand = 2 * pi * rand();

% Calculate the sine-based initialization
% The argument is k*(X + Y) to create a wave diagonal across the grid, 
% plus the random phase offset.
% height_current = A * sin(k * (X + Y) + theta_rand);
height_current = A * randn(simulation_grid_size);


height_time_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);


height_time_series(:, :, 1) = height_current;


E_total = zeros(num_time_steps, 1);
E_bend_history = zeros(num_time_steps, 1);
E_tension_history = zeros(num_time_steps, 1);
E_pin_history = zeros(num_time_steps, 1);
E_area_history = zeros(num_time_steps, 1);


fprintf('Starting simulation with %d steps...\n', num_time_steps);


%% =================  WAVE VECTORS  =================
wavevector = 2*pi/system_length * [0:simulation_grid_size/2, -simulation_grid_size/2+1:-1];
[Qx, Qy] = meshgrid(wavevector, wavevector);
Q2 = Qx.^2 + Qy.^2;
Q2(1,1) = 1e-12; % avoid divide by zero
Q4 = Q2.^2;

% base (diagonal) values for Kq
base_kappa = mean(kappa_map(:));
base_gamma = mean(gamma_map(:));

% de-alias (2/3 rule) mask (centered in FFT ordering)
Nx = simulation_grid_size; Ny = Nx;
cx = floor(Nx/2)+1; cy = floor(Ny/2)+1;
kx_cut = floor((2/3)*(Nx/2));
ky_cut = kx_cut;
mask = zeros(Ny, Nx);
ix = (cx-kx_cut):(cx+kx_cut);
iy = (cy-ky_cut):(cy+ky_cut);
mask(iy, ix) = 1;


% regularized mobility in q-space (Oseen-like), avoid divergence at tiny q
eps_q = (2*pi/system_length)^2 * 1e-8;
Lambda_q = 1 ./ (4 * viscosity * sqrt(Q2 + eps_q));
Lambda_q(1,1) = 0;       % no mobility for q=0 (handle mean separately)


qphys = sqrt(Q2);
qcut = 0.67*(pi/grid_spacing);
filter_q = exp(-(qphys/qcut).^8);   % soft rolloff

figure();
imagesc(qphys);
sprintf('The cutoff is %f',qcut);
figure();
imagesc(filter_q);


% constants for diagnostics or alternative schemes
mobility_det = 1/(viscosity * qcut);  % characteristic (not used in modal update)
% End precompute

%% =================  MAIN LOOP  =================
fig1 = figure('Position', [100, 100, 600, 800], 'Color', 'w');

for step = 2:num_time_steps

    % 1) current Hq (Fourier amplitudes of the field)
    height_current = height_time_series(:,:,step-1);
    Hq = fft2(height_current);   % Nx-by-Nx complex array

    % 2) compute nonlinear deterministic forces in real space (vectorized)
    %    e.g., pinning, external potential, local nonlinear area term if you insist
    %    Recommended: use global-area effective tension instead of local grad_sq*Laplacian
    % real-space fields
    h_x = real(ifft2(1i * Qx .* Hq));
    h_y = real(ifft2(1i * Qy .* Hq));
    Laplacian_h = real(ifft2(-Q2 .* Hq));
    Biharmonic_h = real(ifft2(Q4 .* Hq));

    % Correct bending force with spatially varying kappa:
    tmp1 = kappa_map .* Laplacian_h;               % kappa(x) * ∇^2 h
    %+  real(ifft2( -Q2 .* fft2(tmp2) )); 
    F_bend_full = - real(ifft2( -Q2 .* fft2(tmp1) )) ;

    % Correct tension force with spatially varying gamma:
    gx = gamma_map .* h_x; gy = gamma_map .* h_y;
    F_tens_full = real(ifft2( 1i*(Qx .* fft2(gx) + Qy .* fft2(gy)) ));  % = ∇·(gamma ∇h)

    % Pinning / external (local)
    F_pin = - pinning_potential .* height_current;
    F_ext_local = - V_ext .* height_current;

    % Nonlinear area derivative (if you want the full nonlinear derivative, add here).
    % Otherwise compute global Sigma_area and treat as linear extra tension (goes into Kq)
    grad_sq_map = h_x.^2 + h_y.^2;
    I_grad = sum(grad_sq_map(:)) * (grid_spacing^2);
    %Sigma_area = (area_compressibility / (2 * membrane_area)) * I_grad;
    % full area force (nonlinear piece) — optional, only if you want non-linear correction
    % F_area_nl = -(area_compressibility) * divergence( (|grad h|^2) * grad h )   (compute via spectral/FD if desired)
    % For stability most use Sigma_area * Lap_h as linearized contribution.

    % Uniform pressure: acts on q=0 only; create real-space constant forcing (same everywhere)
    F_press_real = - P_ext * ones(size(height_current));  % careful: this is area force per unit area

    % ------------------- Compose FULL real-space force ------------------------------
    F_full_real = F_bend_full + F_tens_full + F_pin + F_ext_local + F_press_real;
    % Note: Sigma_area * Lap_h is included here only if you want it also to appear as part of F_nl.
    % But we'll include Sigma_area in Kq (preferred), so subtract its linear part below.

    % ------------------- Build linear diagonal Kq using base constants ----------------
    Kq = kappa_map .* Q4 + (gamma_map).* Q2;
    Kq(1,1) = Inf;   % treat q=0 separately

    % ------------------- Extract residual (non-diagonal) forcing ---------------------
    % residual = full force - (linear operator applied in real space)
    % linear_op_real = - base_kappa * ∇^4 h + base_gamma * ∇^2 h + Sigma_area * ∇^2 h + F_press_linear
    % Compute linear_op_real in real space (use same spectral derivatives)
    linear_bend_real = - base_kappa .* Biharmonic_h ; 
    linear_tens_real = base_gamma .* Laplacian_h;
    linear_press_real = - P_ext * ones(size(height_current)); % if you want pressure in linear part
    linear_op_real = linear_bend_real + linear_tens_real + linear_press_real;


    % ------------------- OU modal parameters (use time_step, not step) ----------------
    tau_q = 1 ./ (Lambda_q .* Kq);
    % --- mobility and Kq already computed earlier ---
    % tau_q = 1 ./ (Lambda_q .* Kq);   % you already have this
    % Compute alpha using time_step (NOT step)
    alpha_q = exp(- time_step ./ tau_q);
    alpha_q(isnan(alpha_q)) = 0;

    % Correct discrete variance for MATLAB fft2 normalization (dx = grid_spacing)
    Vq = (boltzmann_const * temperature*mobility_det) ./ (Kq * (grid_spacing^4));
    % clean up numerics
    Vq(isinf(Vq) | isnan(Vq)) = 0;

    % Generate Hermitian Gaussian noise
    Xi = (randn(simulation_grid_size, simulation_grid_size) + 1i*randn(simulation_grid_size, simulation_grid_size));
    Xi = makeHermitian(Xi);   % ensure conjugate symmetry and Nyquist/DC real
    Xi(1,1) = real(Xi(1,1));

    % stochastic increment (per-mode)
    stoch_q = sqrt( Vq .* (1 - alpha_q.^2) ) .* Xi.*filter_q;

    % 11) final modal update
    Hq_new = alpha_q .* Hq  + stoch_q;

    % 12) optional spectral filter and back to real
    Hq_new = Hq_new .* filter_q;
    height_current = real(ifft2(Hq_new));

    % 13) remove mean (control q=0 / volume)
    height_current = height_current - mean(height_current(:));

    % 14) store & visualize
    height_time_series(:,:,step) = height_current;
    if rem(step,5000)==1
        update_visualization(fig1, X_nm, Y_nm,L_nm,height_time_series(:,:,step), time_points(step),step);
    end   
end

%% Generate IRM Images and MoransI Map from the fluctuation map of the simulation


std_dev_time = std(height_time_series,0,3);
figure;
imshow(std_dev_time*1e9,[min(std_dev_time(:)*1e9) max(std_dev_time(:)*1e9)]);
colormap jet;
colorbar;
caxis([min(std_dev_time(:)*1e9) max(std_dev_time(:))*1e9]);
title('Standard Deviation of height(nm)');


ricm_time_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
for step = 1:50:num_time_steps
        height_current = height_time_series(:,:,step);

        % Calculate wave number for the membrane
        wave_number = 2 * pi * refractive_indices(2) / wavelength;

        % Calculate reflection intensities
        reflected_intensity_1 = fc(1, 2, refractive_indices)^2 * base_intensity;
        reflected_intensity_2 = (1 - fc(1, 2, refractive_indices)^2) * base_intensity * fc(2, 3, refractive_indices)^2;

        % Calculate interference terms
        max_intensity = reflected_intensity_1 + reflected_intensity_2 + 2 * sqrt(reflected_intensity_1 * reflected_intensity_2);
        min_intensity = reflected_intensity_1 + reflected_intensity_2 - 2 * sqrt(reflected_intensity_1 * reflected_intensity_2);

        sum_intensity = max_intensity + min_intensity;
        diff_intensity = max_intensity - min_intensity;

        % Calculate reflection coefficients
        reflection_factor = fc(2, 3, refractive_indices) / fc(1, 2, refractive_indices);
        interference_coeff = reflection_factor * (1 - fc(1, 2, refractive_indices)^2);

        % Phase calculations
        phase_shift = 4 * pi * refractive_indices(2) * separation_distance / wavelength;
        reference_height = -(wavelength / (4 * pi * refractive_indices(1))) * ...
            atan((interference_coeff * sin(phase_shift)) / ...
            (1 + interference_coeff * sin(phase_shift)));

        % Final intensity calculation
        intensity_map = sum_intensity/2 - diff_intensity * ...
            cos(2 * wave_number * (height_current - reference_height));
        
        ricm_time_series(:,:,step) = intensity_map;

        % Calculate Moran's I
        kernel_size = 11;
        center = ceil(kernel_size/2);
        [kernel_x, kernel_y] = meshgrid(-(center-1):(kernel_size-center), -(center-1):(kernel_size-center));
        distance_kernel = sqrt(kernel_x.^2 + kernel_y.^2) * grid_spacing;
        epsilon = 1e-10;
        weighting_kernel = 1 ./ (distance_kernel + epsilon);
        weighting_kernel(center, center) = 0;

        global_mean_h = mean(height_current(:));
        mean_centered_h = height_current - global_mean_h;
        weighted_sum_of_neighbors_deviations = conv2(mean_centered_h, weighting_kernel, 'same');
        numerator_map = mean_centered_h .* weighted_sum_of_neighbors_deviations;
        global_variance = var(height_current(:), 1);
        morans_i_map = numerator_map / global_variance;


        % Save intensity frame
        intensity_normalized = mat2gray(intensity_map);
        intensity_uint8 = im2uint8(intensity_normalized);
        filename = sprintf('frame_%04d.png', step);
        imwrite(intensity_uint8, fullfile(output_dir_irm, filename));

        % Save Moran's I frame
        morans_i_normalized = mat2gray(morans_i_map);
        morans_i_uint8 = im2uint8(morans_i_normalized);
        filename = sprintf('MoransI_%04d.png', step);
        imwrite(morans_i_uint8, fullfile(output_dir_moransI, filename));
        
        height_norm = mat2gray(height_time_series(:,:,step));
        height_image = im2uint8(height_norm);
        filename = sprintf('Height_Map_%04d.png',step);
        imwrite(height_image,fullfile(output_dir_curvature_map,filename));
        
        
    
end



[PSD_2D, qx, qy, q_radial, PSD_radial, PSD_qf, freq_Hz, q_selected_list] = compute_membrane_PSDs((height_time_series.*1e9), grid_spacing*1e9, grid_spacing*1e9, time_step);


figure;
imagesc( fftshift(qx(1,:)), fftshift(qy(:,1)), fftshift(log10(abs(PSD_2D))) );
axis xy; xlabel('q_x (rad/m)'); ylabel('q_y (rad/m)');
title('Fourier Space(q_x,q_y)'); colorbar;
% Save figure
saveas(gcf, fullfile(dir, 'Powe Space.png'));
close(gcf);


% =========================================================
% 1) Radially averaged spatial PSD
% =========================================================
figure('Position',[100 100 600 450]);
loglog(q_radial*1e9, PSD_radial, '-o', 'LineWidth',1.5, 'MarkerSize',5);
xlabel('q (1/nm)', 'FontSize', 12);
ylabel('S_{radial}(q) [nm^2]', 'FontSize', 12);
title('Radially Averaged Spatial PSD', 'FontSize', 13);
grid on;
set(gca,'FontSize',11);
xlim([min(q_radial*1e9)*0.8, max(q_radial*1e9)*1.2]);

% Save figure
saveas(gcf, fullfile(dir, 'Radially_Averaged_PSD.png'));
close(gcf);


% =========================================================
% 2) Temporal PSD for selected q-modes
% =========================================================
num_modes = size(PSD_qf,1);
Nfreq     = numel(freq_Hz);
omega     = 2*pi*freq_Hz;     % convert Hz → rad/s

for m = 1:num_modes
    % get the corresponding q-index
    q_idx = q_selected_list(m,:);
    q_value = sqrt(qx(q_idx(1), q_idx(2))^2 + qy(q_idx(1), q_idx(2))^2); % [1/nm]
    q_value_m = q_value * 1e9; % convert [1/nm] → [1/m] for theory
    
    % measured PSD (nm²/Hz)
    PSD_measured = PSD_qf(m,:);


    % ---- Plotting ----
    figure('Position',[100 100 600 450]);
    loglog(freq_Hz, PSD_measured, 'b-', 'LineWidth',1.5, ...
        'DisplayName', 'Measured PSD');

    xlabel('Frequency f (Hz)', 'FontSize', 12);
    ylabel('PSD (nm^2/Hz)', 'FontSize', 12);
    title(sprintf('Temporal PSD for Mode |q| = %.2e [1/m]', q_value_m), 'FontSize', 13);
    legend('Location','southwest');
    grid on;
    set(gca,'FontSize',11);
    xlim([freq_Hz(2), max(freq_Hz)/2]);

    % Save each plot
    filename = sprintf('Temporal_PSD_q%.2e.png', q_value_m);
    saveas(gcf, fullfile(dir, filename));
    close(gcf);
end


function [PSD_2D, qx, qy, q_radial, PSD_radial, PSD_qf, freq_Hz, q_selected_list] = compute_membrane_PSDs(height_series_nm, dx_nm, dy_nm, dt_s)
% =========================================================================
% COMPUTE_MEMBRANE_PSDS
% -------------------------------------------------------------------------
% Computes 2D spatial and spatio-temporal Power Spectral Densities (PSDs)
% for a fluctuating membrane height field h(x,y,t).
%
% INPUTS:
%   height_series_nm : 3D array [Ny x Nx x Nt]  - height field in nm
%   dx_nm, dy_nm     : pixel sizes in nm
%   dt_s             : time step between frames in seconds
%
% OUTPUTS:
%   PSD_2D           : 2D spatial PSD averaged over time [nm²]
%   qx, qy           : 2D wavevector grids [1/nm]
%   q_radial         : radial wavenumber bins [1/nm]
%   PSD_radial       : radially averaged PSD [nm²]
%   PSD_qf           : spatio-temporal PSD for selected q values [nm²/Hz]
%   freq_Hz          : temporal frequency array [Hz]
%   q_selected_list  : list of selected (y,x) indices used for temporal PSDs
%
% UNITS:
%   PSD_2D(qx,qy): nm²
%   PSD_qf(q,f):   nm²/Hz
%
% PHYSICAL NORMALIZATION:
%   sum(PSD_2D)*ΔqxΔqy ≈ variance(h)
%
% =========================================================================

[Ny, Nx, Nt] = size(height_series_nm);

% -------------------------------------------------------------------------
% 1. Geometry & Wavevector grids
% -------------------------------------------------------------------------
Lx_nm = Nx * dx_nm;
Ly_nm = Ny * dy_nm;

kx = fftshift( (2*pi/Lx_nm) * (-Nx/2:Nx/2-1) );  % [1/nm]
ky = fftshift( (2*pi/Ly_nm) * (-Ny/2:Ny/2-1) );
[qx, qy] = meshgrid(kx, ky);
q_magnitude = sqrt(qx.^2 + qy.^2);

% -------------------------------------------------------------------------
% 2. Accumulators
% -------------------------------------------------------------------------
PSD_sum = zeros(Ny, Nx);
q_selected_list = select_representative_q_indices(Nx, Ny,10,2);
num_q_selected = size(q_selected_list,1);
fourier_mode_time_series = zeros(num_q_selected, Nt);  % complex

% -------------------------------------------------------------------------
% 3. Loop over time frames
% -------------------------------------------------------------------------
for t = 1:Nt
    h_nm = height_series_nm(:,:,t);
    F_nm = fftshift(fft2(h_nm)) * dx_nm * dy_nm;  % physical scaling
    PSD_sum = PSD_sum + abs(F_nm).^2;

    % Store selected q-modes
    for q_idx = 1:num_q_selected
        iy = q_selected_list(q_idx,1);
        ix = q_selected_list(q_idx,2);
        fourier_mode_time_series(q_idx, t) = F_nm(iy, ix);
    end
end

% -------------------------------------------------------------------------
% 4. Time-averaged 2D PSD
% -------------------------------------------------------------------------
area_nm2 = Lx_nm * Ly_nm;
PSD_2D = PSD_sum / (area_nm2 * Nt);  % [nm²]

% -------------------------------------------------------------------------
% 5. Radial averaging
% -------------------------------------------------------------------------
num_bins = 100;
q_edges = linspace(0, max(q_magnitude(:)), num_bins + 1);
q_centers = 0.5 * (q_edges(1:end-1) + q_edges(2:end));
PSD_radial = zeros(1, num_bins);
bin_counts = zeros(1, num_bins);

for iy = 1:Ny
    for ix = 1:Nx
        q_val = q_magnitude(iy,ix);
        bin_idx = find(q_val >= q_edges(1:end-1) & q_val < q_edges(2:end), 1);
        if isempty(bin_idx), continue; end
        PSD_radial(bin_idx) = PSD_radial(bin_idx) + PSD_2D(iy,ix);
        bin_counts(bin_idx) = bin_counts(bin_idx) + 1;
    end
end

valid = bin_counts > 0;
PSD_radial(valid) = PSD_radial(valid) ./ bin_counts(valid);
q_radial = q_centers;

% -------------------------------------------------------------------------
% 6. Temporal PSD for selected q-values
% -------------------------------------------------------------------------
nfft = 2^nextpow2(Nt);
freq_Hz = (0:nfft/2-1) / (nfft*dt_s);
PSD_qf = zeros(num_q_selected, length(freq_Hz));

for q_idx = 1:num_q_selected
    signal = squeeze(fourier_mode_time_series(q_idx,:));
    signal = signal - mean(signal);
    window = hann(Nt)';
    sig_win = signal .* window;
    spec = fft(sig_win, nfft);
    Pxx = (2*dt_s / (Nt * sum(window.^2))) * abs(spec(1:nfft/2)).^2;
    PSD_qf(q_idx,:) = Pxx;  % [nm²/Hz]
end

end

% -------------------------------------------------------------------------
% Helper: choose representative q indices (low q  and isotropic)
% -------------------------------------------------------------------------
function q_indices = select_representative_q_indices(Nx, Ny, max_radius, mesh_step)
    % ---------------------------------------------------------------------
    % Select low-q (long wavelength) Fourier modes near center
    % using a regular mesh with configurable spacing.
    %
    % Inputs:
    %   Nx, Ny      - image dimensions
    %   max_radius  - radius (in pixels) from DC to include
    %   mesh_step   - pixel spacing in Fourier grid (e.g. 2 = every 2 pixels)
    %
    % Output:
    %   q_indices   - [Nq x 2] array of (y,x) indices
    % ---------------------------------------------------------------------

    if nargin < 3
        max_radius = 10;  % default radius (pixels)
    end
    if nargin < 4
        mesh_step = 2;    % default step of 2 pixels
    end

    cx = floor(Nx/2) + 1;
    cy = floor(Ny/2) + 1;

    [X, Y] = meshgrid(1:mesh_step:Nx, 1:mesh_step:Ny);
    R = sqrt((X - cx).^2 + (Y - cy).^2);

    % Mask for low-q region, excluding DC (optional)
    mask = (R > 0) & (R <= max_radius);

    iy = Y(mask);
    ix = X(mask);

    q_indices = [iy(:), ix(:)];

    % Ensure unique + in-bound indices
    q_indices(:,1) = mod(q_indices(:,1)-1, Ny) + 1;
    q_indices(:,2) = mod(q_indices(:,2)-1, Nx) + 1;
    q_indices = unique(q_indices, 'rows', 'stable');
end




%% ================= VISUALIZATION FUNCTION =================
function  update_visualization (fig, X_nm, Y_nm, L_nm, h_total, current_time,step)
    figure(fig);
    clf;

    h_nm = h_total * 1e9;

    % Main surface plot
    subplot(2, 2, [1, 2]);
    % Use the full 2D arrays for X, Y, and Z
    surf(X_nm, Y_nm, h_nm, 'EdgeColor', 'none');
    title(sprintf('Membrane Height at t = %.2f ms', current_time*1e3));
    xlabel('x (nm)');
    ylabel('y (nm)');
    zlabel('Height (nm)');
    colormap(turbo);
    colorbar;
    view(-30, 60);
    

    % Top view
    subplot(2,2,3);
    contourf(X_nm, Y_nm, h_nm, 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Height contour map');
    
    % Height distribution
    subplot(2,2,4);
    histogram(h_nm(:), 50, 'Normalization', 'pdf');
    title('Height Distribution');
    xlabel('Height (nm)'); ylabel('Probability Density');
    grid on;
    drawnow;
    
 
end


function frensel_coefficient = fc(index1,index2,refractive_index)
n_1 = refractive_index(index1);
n_2 = refractive_index(index2);
frensel_coefficient=(n_1-n_2)/(n_1+n_2);
end

function X = makeHermitian(X)
    [Ny, Nx] = size(X);

    % Conjugate symmetry
    for y = 2:Ny
        for x = 2:Nx
            X(Ny - y + 2, Nx - x + 2) = conj(X(y, x));
        end
    end

    % Force Nyquist frequencies and DC to be real
    X(1,1)           = real(X(1,1));          % DC
    X(1,Nx/2+1)      = real(X(1,Nx/2+1));     % Nyquist x-axis
    X(Ny/2+1,1)      = real(X(Ny/2+1,1));     % Nyquist y-axis
    X(Ny/2+1,Nx/2+1) = real(X(Ny/2+1,Nx/2+1));% Nyquist corner
end

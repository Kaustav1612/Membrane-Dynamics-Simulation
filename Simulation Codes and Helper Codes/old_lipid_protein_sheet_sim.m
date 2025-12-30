% Membrane Dynamics Simulation with Physical Mode Evolution and Advanced Visualization
clear; close all; clc;

%% ================= PHYSICAL PARAMETERS =================
simulation_grid_size = 60;       % Reduced for faster testing
simulated_area =1e-12;           % 5μm x 5μm membrane (15 μm²)
membrane_area = 1e-12;            % Membrane area should match simulated area
viscosity = 140;                 % Water viscosity (Pa·s)
temperature = 310;                % Temperature (K)
boltzmann_const = 1.38e-23;       % Boltzmann constant (J/K)
bending_rigidity = 18 * boltzmann_const * temperature; % ~20 kT
time_step = 10e-6;                 % MUCH smaller time step (20 ms)
total_simulation_time = 5;     % 10 ms total timefTpo
surface_tension = 500e-6;           % MUCH lower tension (500 pN/um)                   % Shear Modulus (N/m)                  
old_dir = pwd;
dir = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein';
output_dir_irm = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein\results_irm';            
output_dir_moransI = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein\results_moransI';
output_dir_bending_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein\bending_map';
output_dir_lipid_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein\results_lipid_map';
output_dir_protein_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein\results_protein_map';
output_dir_curvature_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\raft_loving_protein\results_curvature_map';
pinning_strength = 30e-6;          % Moderate pinning (N/m²)


refractive_indices = [1.33, 1.48, 1.35];  % [air,media,membrane,cytosol]
separation_distance = 5e-9;
base_intensity = 0.1;
wavelength = 546e-9;
V_ext = ones(simulation_grid_size)*5e-6;  % Example: external potential map
P_ext = ones(simulation_grid_size)*0.1e-6;                         % External pressure term
low_kappa = 15;
high_kappa = 20;
low_tension = 250e-6;
high_tension  = 750e-6;
beta_phi   = 0.5;              % strong stiffening by proteins
beta_psi   = 0.5;              % moderate stiffening by rafts
C_p        = 1/(1e-6);             %(protein curvature)
C_l        = 1/(2e-6);             % (lipid curvature)

% Protein clustering properties
N = simulation_grid_size;
A_hat = 20;     % interaction parameter
S_hat = 200;     % saturation denstiy parameter
phi0_mean = 0.2;    % mean area fraction
sigma_init = 0.005;  % noise strength
phi = phi0_mean + imgaussfilt(sigma_init*randn(N,N),[5,5]);
phi = min(max(phi,1e-6),1-1e-6);   % enforce bounds
phi_eq = phi0_mean + 0.02*randn(size(phi));
phi_eq = max(min(phi_eq, 0.01), 0.05); % keep within [0.01, 0.1]

% Lipid heterogenity properties 
N = simulation_grid_size;
alpha = 0.1;    
J = 0.85;
lambda = 0.01;
psi0_mean = 0.2;    % mean area fraction
sigma_init = 0.05;  % noise strength
psi = psi0_mean + imgaussfilt(sigma_init*randn(N,N),[5,5]);
psi = min(max(psi,1e-6),1-1e-6);   % enforce bounds


C0_map = C_p.*phi + C_l.*psi;

mobility_phi = 0.1e-6;   % Protein mobility (fast)
mobility_psi = 0.1e-8;     % Lipid mobility (slow)

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

% Variation amplitude (±10% variation)
kappa_variation = 0.05 * base_kappa;
gamma_variation = 0.05 * base_gamma;


rng(42); % reproducibility

% white noise
noise_kappa = randn(simulation_grid_size, simulation_grid_size);
noise_gamma = randn(simulation_grid_size, simulation_grid_size);

% smooth to introduce correlation
corr_length = 2; % adjust to control patch size
noise_kappa_smooth = imgaussfilt(noise_kappa, corr_length);
noise_gamma_smooth = imgaussfilt(noise_gamma, corr_length);

% normalise to [-1,1]
noise_kappa_smooth = noise_kappa_smooth ./ max(abs(noise_kappa_smooth(:)));
noise_gamma_smooth = noise_gamma_smooth ./ max(abs(noise_gamma_smooth(:)));

% apply variation (now can go above or below base)
kappa_map = base_kappa + kappa_variation * noise_kappa_smooth;
gamma_map = base_gamma + gamma_variation * noise_gamma_smooth;

kappa_map = imgaussfilt(kappa_map, 0.8);  % sigma in pixels; adjust
gamma_map = imgaussfilt(gamma_map, 0.8);

num_time_steps = round(total_simulation_time / time_step);


%% ================= VISUALIZATION =================
figure();
% Plot kappa map
imagesc(1:simulation_grid_size, 1:simulation_grid_size, kappa_map/(boltzmann_const*temperature));
colormap jet;
axis square; colorbar; 
caxis([low_kappa high_kappa]);
title('Bending Rigidity (κ) Map in kBT');
xlabel('X (pixels)'); ylabel('Y (pixels)');



% Plot gamma map
figure();
imagesc(1:simulation_grid_size, 1:simulation_grid_size, gamma_map*1e6);
colormap jet;
caxis([low_tension*1e6 high_tension*1e6]);
axis square; colorbar;
title('Surface Tension (γ) Map in pN/um');
xlabel('X (pixels)'); ylabel('Y (pixels)');




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

time_points = (0:num_time_steps-1) * time_step;

% Initialize arrays
height_current = 1e-10 * randn(simulation_grid_size);

height_time_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);

height_time_series(:, :, 1) = height_current;
fprintf('Starting simulation with %d steps...\n', num_time_steps);

% phi initialization (example)
protein_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
lipid_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
curvature_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
kappa_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);



% Heterogentity Map
protein_map_series(:,:,1) = phi;
lipid_map_series (:,:,1) = psi;
curvature_map_series (:,:,1) = C0_map;


%% =================  WAVE VECTORS  =================
wavevector_height = 2*pi/system_length * [0:simulation_grid_size/2, -simulation_grid_size/2+1:-1];
[Qx, Qy] = meshgrid(wavevector_height, wavevector_height);
Q2 = Qx.^2 + Qy.^2;
Q2(1,1) = 1e-12;
Q4 = Q2.^2;

wavevector_protein = 2*pi/(system_length*1e6) * [0:simulation_grid_size/2, -simulation_grid_size/2+1:-1];
[Qx_p, Qy_p] = meshgrid(wavevector_protein, wavevector_protein);
Q2_p = Qx_p.^2 + Qy_p.^2;     % Laplacian operator in Fourier space
Q2_p(1,1) = 1e-12;  % Explicit DC component to be very low
Q4_p = Q2_p.^2;             % Biharmonic if needed

wavevector_lipid = 2*pi/(system_length*1e6) * [0:simulation_grid_size/2, -simulation_grid_size/2+1:-1];
[Qx_l, Qy_l] = meshgrid(wavevector_lipid, wavevector_lipid);
Q2_l = Qx_l.^2 + Qy_l.^2;     % Laplacian operator in Fourier space
Q2_l(1,1) = 1e-12; % Explicit DC component to be very low
Q4_l = Q2_l.^2;             % Biharmonic if needed



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
fig2 = figure('Position', [100, 100, 800, 800], 'Color', 'w');
for step = 2:num_time_steps

    % 1) current Hq (Fourier amplitudes of the field)
    height_current = height_time_series(:,:,step-1);
    Hq = fft2(height_current);   % Nx-by-Nx complex array

    % === Gradients and Laplacians ===
    
    % === Gradiens and Laplacians for heights === %
    h_x = real(ifft2(1i * Qx .* Hq));
    h_y = real(ifft2(1i * Qy .* Hq));
    Laplacian_h = real(ifft2(-Q2 .* Hq));
    Biharmonic_h = real(ifft2(Q4 .* Hq));
    
    % === Gradients and Laplacians for proteins ===%
    phi = protein_map_series(:,:,step-1);
    phi_q = fft2(phi);
    grad_phi_x = real(ifft2(1i * Qx_p .* phi_q));
    grad_phi_y = real(ifft2(1i * Qy_p .* phi_q));
    grad_sq_phi = grad_phi_x.^2 + grad_phi_y.^2;    
    Laplacian_phi = real(ifft2(-Q2_p .* phi_q));
    Biharmonic_phi = real(ifft2(Q4_p .* phi_q));
    
    % === Gradients and Laplacians for lipids ===%
    psi = lipid_map_series(:,:,step-1);
    psi_q = fft2(psi);
    grad_psi_x = real(ifft2(1i * Qx_p .* psi_q));
    grad_psi_y = real(ifft2(1i * Qy_p .* psi_q));
    grad_sq_psi = grad_psi_x.^2 + grad_psi_y.^2;
    Laplacian_psi = real(ifft2(-Q2_l .* psi_q));
    Biharmonic_psi = real(ifft2(-Q4_l.* psi_q));

   

    % Cahn Hillard Type of Protein Diffusion and Dynamics

    coeff1 = 1 ./ (1 - phi) - A_hat .* phi;
    coeff2 = 1 ./ (1 - phi).^2 - A_hat;
    coeff3 = (A_hat ./ (2*S_hat));
    term1 = Laplacian_phi .* coeff1;
    term3 = -phi .* (coeff3 .* Biharmonic_phi);
    term2 = grad_sq_phi .* coeff2;
    
    % Raft Protein Coupling J
    % Cytoskeleton unbinding time Tau
    % M_phi, M_psi are mobilities (lipid and protein respectively)
    
    % ---------------- Chemical potentials ----------------
    % Lipid order parameter (raft vs non-raft)
    mu_psi = alpha .* (psi.^3 - psi) - lambda .* Laplacian_psi - J*(psi-phi);
    

    % Protein area fraction (Cahn–Hilliard(from Rangamani Paper)+ nocoupling + height term + no activity)
    mu_phi = term1 + term3 + term2 ...
        - kappa_map.*(Laplacian_h - C_p.*phi);


    % ---------------- Time updates ----------------
    % Protein evolution (Cahn–Hilliard type)
    phi_new = phi + time_step * mobility_phi .* real(ifft2(-Q2_p .* fft2(mu_phi)));

    % Lipid evolution (Cahn–Hilliard)
    psi_new = psi + time_step * mobility_psi .* real(ifft2(-Q2_l .* fft2(mu_psi)));
   
    % Bound enforcement (since phi and psi are bounded order parameters)
    phi_new = min(max(phi_new, 0), 1);
    psi_new = min(max(psi_new, 0), 1);
    

    % Store updated fields
    protein_map_series(:,:,step) = phi_new;
    lipid_map_series(:,:,step)   = psi_new;
    

     % --- Update kappa map (spatially and temporally varying) ---
    kappa_map = base_kappa .* (1 + beta_phi * phi_new + beta_psi * psi_new);
    kappa_map_series(:,:,step) = kappa_map;
    
    protein_map_series(:,:,step) = phi_new;
    % Spontaneous curvature coupling
    C0_map = C_p * phi_new + C_l* psi_new;
   
    % Correct bending force with spatially varying kappa:
    field = kappa_map .* (Laplacian_h - C0_map);
    F_bend_full = real(ifft2(-Q2 .* fft2(field)));  

    % Correct tension force with spatially varying gamma:
    gx = gamma_map .* h_x;
    gy = gamma_map .* h_y;
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
    kappa_mean = mean(kappa_map(:));
    gamma_mean = mean(gamma_map(:));
    
    Kq = kappa_mean .* Q4 + gamma_mean .* Q2;
    Kq(1,1) = Inf; % avoid division by zero for q=0


    % ------------------- Extract residual (non-diagonal) forcing ---------------------
    % residual = full force - (linear operator applied in real space)
    % linear_op_real = - base_kappa * ∇^4 h + base_gamma * ∇^2 h + Sigma_area * ∇^2 h + F_press_linear
    % Compute linear_op_real in real space (use same spectral derivatives)
    linear_bend_real = - base_kappa .* Biharmonic_h ; 
    linear_tens_real = base_gamma .* Laplacian_h;
    linear_press_real = - P_ext * ones(size(height_current)); % if you want pressure in linear part
    linear_op_real = linear_bend_real + linear_tens_real + linear_press_real;

    % Residual forcing (what goes into F_nl_real)
    F_nl_real_residual = F_full_real - linear_op_real;
    % zero out tiny numerical noise
    F_nl_real_residual(abs(F_nl_real_residual) < 1e-20) = 0;

    % ------------------- Transform residual into q-space (dealias / mask optionally) -
    Fq_nl = fft2(F_nl_real_residual);
    Fq_nl = Fq_nl .* filter_q;   % if you use a de-alias mask (recommended)

    % ------------------- OU modal parameters (use time_step, not step) ----------------
    tau_q = 1 ./ (Lambda_q .* Kq);
    % --- mobility and Kq already computed earlier ---
    % tau_q = 1 ./ (Lambda_q .* Kq);   % you already have this
    % Compute alpha using time_step (NOT step)
    alpha_q = exp(- time_step ./ tau_q);
    alpha_q(isnan(alpha_q)) = 0;

    % Correct discrete variance for MATLAB fft2 normalization (dx = grid_spacing)
    Vq = mobility_det* (boltzmann_const * temperature) ./ (Kq * (grid_spacing^4));
    % clean up numerics
    Vq(isinf(Vq) | isnan(Vq)) = 0;

    % 9) deterministic modal forcing from residual
    det_nl_q = (1 - alpha_q) .* (Fq_nl ./ Kq);
    det_nl_q(isnan(det_nl_q) | isinf(det_nl_q)) = 0;

    % Generate Hermitian Gaussian noise
    Xi = (randn(simulation_grid_size, simulation_grid_size) + 1i*randn(simulation_grid_size, simulation_grid_size));
    Xi = makeHermitian(Xi);   % ensure conjugate symmetry and Nyquist/DC real
    Xi(1,1) = real(Xi(1,1));

    % stochastic increment (per-mode)
    stoch_q = sqrt( Vq .* (1 - alpha_q.^2) ) .* Xi.*filter_q;

    % 11) final modal update
    Hq_new = alpha_q .* Hq + det_nl_q + stoch_q;

    % 12) optional spectral filter and back to real
    Hq_new = Hq_new .* filter_q;
    height_current = real(ifft2(Hq_new));

    % 13) remove mean (control q=0 / volume)
    height_current = height_current - mean(height_current(:));

    % 14) store & visualize
    height_time_series(:,:,step) = height_current;
  
    % Update visualization less frequently
    if mod(step,5000) == 1
        update_visualization(fig1,fig2, X_nm, Y_nm,L_nm,height_current, time_points(step),lipid_map_series,kappa_map_series,protein_map_series,step,C0_map);
        drawnow;
    end
end


%%
cd(dir);
% Setup for video
video_filename = 'membrane_dynamics_movie(pinned).mp4';
video1 = VideoWriter(video_filename, 'MPEG-4');
video1.FrameRate = 100;
open(video1);
figure();  % Create and reuse figure window
h_range = max(abs(height_time_series(:)));

cmin = -h_range*1e9;
cmax=  h_range*1e9 ;
for t = 1:100:num_time_steps
    clf;
    h = height_time_series(:,:,t)*1e9;
    surf(X_nm, Y_nm, h, 'EdgeColor','none');
    hold on;

    % Use actual time value in microseconds
    title(sprintf('Membrane Height (t = %.2f ms)', time_points(t)*1e3));
    xlabel('x/L'); ylabel('y'); zlabel('h (nm)');
    view(-30, 60);
    colormap(turbo);
    hcb = colorbar('Location','eastoutside');
    hcb.Label.String = 'Height (nm)';
    caxis([cmin cmax]);
    zlim([cmin cmax]);
    light('Position', [1 1 1], 'Style', 'infinite');

    drawnow;
    frame = getframe(gcf);

    writeVideo(video1, frame);
end

close(video1);
cd(old_dir);

%% Generate IRM Images and MoransI Map from the fluctuation map of the simulation

std_dev_time = std(height_time_series,0,3);
figure;
imshow(std_dev_time*1e9,[min(std_dev_time(:)*1e9) max(std_dev_time(:)*1e9)]);
colormap jet;
colorbar;
caxis([min(std_dev_time(:)*1e9) max(std_dev_time(:))*1e9]);
title('Standard Deviation of height(nm)');
filename = sprintf('SD_Time.png');
full_path = fullfile(dir, filename);
saveas(gcf, full_path);
close;

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
        
        phi_norm = mat2gray(protein_map_series(:,:,step));
        phi_image = im2uint8(phi_norm);
        filename = sprintf('Protein_Density_Map_%04d.png',step);
        imwrite(phi_image,fullfile(output_dir_protein_map,filename));
        
        
        a = lipid_map_series(:,:,step);
        binary_lipid_map_series(:,:,step) = a > 0.01;
        psi_norm = mat2gray(lipid_map_series(:,:,step));
        phi_image = im2uint8(psi_norm);
        filename = sprintf('Lipid_Density_Map_%04d.png',step);
        imwrite(phi_image,fullfile(output_dir_lipid_map,filename));
        
        
        height_norm = mat2gray(height_time_series(:,:,step)*1e9);
        height_image = im2uint8(height_norm);
        filename = sprintf('Height_Map_%04d.png',step);
        imwrite(height_image,fullfile(output_dir_curvature_map,filename));

end

%% Finding the PSD of the simulation data

[S2D, qx, qy, q_r, S_radial, S_qw, omega, q_list] = compute_PSDs(height_time_series, grid_spacing, grid_spacing, time_step);

% 1) 2D PSD heatmap (log scale)
figure;
imagesc( fftshift(qx(1,:)), fftshift(qy(:,1)), fftshift(log10(abs(S2D))) );
axis xy; xlabel('q_x (rad/m)'); ylabel('q_y (rad/m)');
title('log10 S(q_x,q_y)'); colorbar;

% 2) radial PSD
figure;
loglog(q_r, S_radial, '-o');
xlabel('q (rad/m)'); ylabel('S_{radial}(q)');
title('Circularly averaged PSD');
filename = sprintf('Circularly averaged PSD.png');
full_path = fullfile(dir, filename);
saveas(gcf, full_path);
close;

% 3) temporal PSD for selected q-modes
for qq = 1:size(S_qw,1)
    figure;
    loglog(omega, S_qw(qq,:));
    xlim([omega(2), max(omega)/2]);
    title(sprintf('Mode %d temporal PSD', qq));
    xlabel('\omega (rad/s)'); ylabel('S(\omega)');
    filename = sprintf('Temporal_PSD_mode%d.png',qq);
    full_path = fullfile(dir, filename);
    saveas(gcf, full_path);
    close;
end



function [S2D, qx, qy, q_r, S_radial, S_qw, omega, q_list] = compute_PSDs(h_ts, dx, dy, dt)
    [Nx, Ny, Nt] = size(h_ts);
    % Note: in your code Nx==Ny==simulation_grid_size; but handle general Nx,Ny.
    Nx = size(h_ts,2);
    Ny = size(h_ts,1); % if orientation differs, ensure consistent: here assume (Ny x Nx x Nt)
    [Ny, Nx, Nt] = size(h_ts);

    Lx = Nx * dx;
    Ly = Ny * dy;
    A = Lx*Ly;
    
    % wavevector grids (MATLAB fft ordering)
    kx = (2*pi/Lx) * [0:(Nx/2-1) (-Nx/2):-1];
    ky = (2*pi/Ly) * [0:(Ny/2-1) (-Ny/2):-1];
    [Qx, Qy] = meshgrid(kx, ky);
    q2 = Qx.^2 + Qy.^2;
    q_mag = sqrt(q2);

    % accumulate |H(q)|^2 over time
    H2_acc = zeros(Ny, Nx);
    % store Hq(t) for spatio-temporal PSD for a subset of q (avoid huge memory)
    % choose a few representative q vectors (e.g., along kx axis and some radial rings)
    q_indices = select_q_indices(Nx, Ny); % helper to choose list of (iy,ix) pairs
    n_qsel = size(q_indices,1);
    Hq_time = zeros(n_qsel, Nt);  % complex

    for tt = 1:Nt
        h = h_ts(:,:,tt);
        H = fft2(h);                     % MATLAB fft2: no forward normalization
        H2_acc = H2_acc + abs(H).^2;     % accumulate power per mode

        % collect chosen mode coefficients over time
        for qq = 1:n_qsel
            iy = q_indices(qq,1); ix = q_indices(qq,2);
            Hq_time(qq, tt) = H(iy, ix);
        end
    end

    % Time-averaged 2D PSD (normalization explained below)
    % Using normalization: S2D(q) = (dx*dy)/(Nx*Ny*Nt) * sum_t |H(q,t)|^2
    S2D = (dx*dy) / (Nx * Ny * Nt) * H2_acc;

    % Circular / radial average
    qmax = max(q_mag(:));
    nbins = floor(min(Nx,Ny)/2);
    q_edges = linspace(0, qmax, nbins+1);
    q_centers = 0.5*(q_edges(1:end-1) + q_edges(2:end));
    S_radial = zeros(1, nbins);
    counts = zeros(1, nbins);

    for iy = 1:Ny
        for ix = 1:Nx
            qv = q_mag(iy,ix);
            % find bin
            bin = find(qv >= q_edges(1:end-1) & qv < q_edges(2:end), 1);
            if isempty(bin)
                continue
            end
            S_radial(bin) = S_radial(bin) + S2D(iy,ix);
            counts(bin) = counts(bin) + 1;
        end
    end
    nonzero = counts>0;
    S_radial(nonzero) = S_radial(nonzero) ./ counts(nonzero);
    q_r = q_centers;

    % ---- Spatio-temporal PSD for selected q's: do temporal FFT of Hq_time ----
    % Option: detrend or window time series before FFT to reduce leakage
    nfft = 2^nextpow2(Nt);
    df = 1/(dt * nfft);
    freq = (0:(nfft-1)) * df;                 % Hz
    omega = 2*pi*freq;                        % rad/s
    S_qw = zeros(n_qsel, nfft);               % power spectral density per mode vs freq

    for qq = 1:n_qsel
        hq_t = squeeze(Hq_time(qq,:));        % complex time series (units inherited from fft2)
        % optional: remove mean
        hq_t = hq_t - mean(hq_t);
        % window
        w = hann(Nt)';
        hq_tw = hq_t .* w;
        H_t = fft(hq_tw, nfft);
        Pxx = (dt / (Nt * sum(w.^2))) * abs(H_t).^2;  % one-sided not done; keep full-spectrum
        S_qw(qq, :) = Pxx;
    end

    % return grid and list of selected q indices for reference
    qx = Qx; qy = Qy;
    q_list = q_indices;
end

%% helper: choose representative q indices (center + few rings + kx axis)
function q_indices = select_q_indices(Nx, Ny)
    cx = floor(Nx/2)+1; cy = floor(Ny/2)+1; % MATLAB 1-based center
    % pick small set: center, few along +kx and +ky and diagonals
    picks = [cy, cx; cy, cx+1; cy, cx+2; cy+1, cx; cy+2, cx; cy+1, cx+1; cy+2, cx+2];
    % clip indices within [1,N]
    picks(:,1) = mod(picks(:,1)-1, Ny)+1;
    picks(:,2) = mod(picks(:,2)-1, Nx)+1;
    q_indices = unique(picks, 'rows', 'stable');
end





%% ================= VISUALIZATION FUNCTION =================
function update_visualization(fig1,fig2,X_nm, Y_nm, L_nm, h_total, current_time,lipid_map_series,kappa_map_series,protein_map_series,step,C0_map)
    figure(fig1);
    clf;
    
    % Assuming X_nm and Y_nm are your 2D meshgrid arrays, and L_nm is the domain length.
    % Convert height to nm
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
   

    figure(fig2);
    clf;
     % Cross-section
    subplot(3,1,1);
    contourf(X_nm, Y_nm, protein_map_series(:,:,step), 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Protein map');


    subplot(3,1,2);
    contourf(X_nm, Y_nm, kappa_map_series(:,:,step), 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Bending rigidity map');
    
    a = lipid_map_series(:,:,step);
    binary_lipid_map_series(:,:,step) = a > 0.015;
    subplot(3,1,3);
    contourf(X_nm, Y_nm, C0_map, 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Curvature map');
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

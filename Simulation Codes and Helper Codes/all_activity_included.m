% Membrane Dynamics Simulation with Physical Mode Evolution and Advanced Visualization
clear; close all; clc;

%% ================= PHYSICAL PARAMETERS =================
simulation_grid_size = 60;       % Reduced for faster testing
simulated_area =5e-12;        % 5μm x 5μm membrane (15 μm²)
membrane_area = 5e-12;            % Membrane area should match simulated area
viscosity = 10;                 % Water viscosity (Pa·s)
temperature = 310;                % Temperature (K)
boltzmann_const = 1.38e-23;       % Boltzmann constant (J/K)
bending_rigidity = 12 * boltzmann_const * temperature; % ~20 kT
time_step = 1e-3;                 % MUCH smaller time step (20 ms)
total_simulation_time = 100;     % 200 s total time
surface_tension = 40e-6;           % high tension (60 pN/um)
% shear_modulus = ;
% confinement_matrix = ;

old_dir = pwd;
dir = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\active\pinned\protrusion';
% output_dir_bending_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\low_curvature\results_bending_map';
% output_dir_tension_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\low_curvature\results_tension_map';
% output_dir_lipid_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\low_curvature\results_lipid_map';
% output_dir_protein_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\low_curvature\results_protein_map';
% output_dir_curvature_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\low_curvature\results_curvature_map';
% output_dir_cyto_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\passive\low_curvature\results_cyto_map';
pinning_strength = 5000;          % Confinement Parameter
low_kappa = 12;
high_kappa = 20;
low_tension = 35e-6;
high_tension  = 75e-6;
K_A = 0;


% --- Protein and lipid coupling parameters ---

C_p = 0.009e9;                % protein curvature [1/m]
C_ld = 0.0008e9;               % lipid curvature [1/m]
C_lo = 0.0005e9;

beta_phi = 0.8;              % protein stiffening strength
beta_psi = 0.2;              % lipid stiffening strength
gamma_phi =  -0.5;
gamma_psi  = 0.5;

rate_exo_per_m2 = 0.005;
rate_endo_per_m2 = 0.01;

phi_eq = 0.35;                  % equilibrium protein order
psi_eq = 0;                  % equilibrium lipid order
tau_phi = 0.5;                 % remodeling timescale for phi [s]
tau_psi = 0.25;               % relaxation timescale for psi [s]


N = simulation_grid_size;
A_hat = 40;     % interaction parameter
S_hat = 400;     % saturation denstiy parameter
alpha = 5;    
J = 0.75;       % J increases attachment in rafts, 
J_l = 0.0001; 
lambda = 0.01;

% --- Protein field phi ---
phi = 0.5 * randn(N, N);  % small random fluctuations around 0
phi = imgaussfilt(phi, 1);   % spatially smooth
phi(phi > 1) = 1; phi(phi < -1) = -1;

% --- Lipid order parameter psi ---
psi = 0.05 * randn(N, N);  % small fluctuations
psi = imgaussfilt(psi, 1);
psi(psi > 1) = 1; psi(psi < -1) = 0;

mobility_phi = 0.1e-12;   % Protein mobility (fast)
mobility_psi = 1e-12;     % Lipid mobility (slow)

% Cytoskeleton turnover parameters
tau_cs = 300;       % cytoskeleton remodelling time  in s
s_activity = 0.3;     % activity factor (0..1). Larger s => more assembly relative to detachment
P_max = 1e4; 
tau_withdrawal = 500; % cytoskeleton force correlation time in s

% Protein production/decay kinetics (protein density model)
k_on_base = 0.01;      % baseline attachment/production rate (per time unit)
k_off = 0.0075;         % baseline detachment/decay rate (per time unit)
max_local_density = 0.5; % saturating density (arbitrary units) to prevent runaway growth


% Gate/threshold for cytoskeleton enabling clusters (optional)
cs_threshold = 0.5;   % if you want a binary gating on cs_presence probability map

% Mobility gating (so proteins are mobile/cluster only where cytoskel present)
mobility_phi_max = mobility_phi; % store original mobility
% Intended logic:
C_l = C_ld * ones(size(psi));
C_l(psi > 0.5) = C_lo;      
C0_map = C_p.*phi + C_l;

tau_a = 1000*time_step;        % Persistence time of activity
D_a = 1*boltzmann_const*temperature;          % Strength of active noise 

sprintf('Activity of the system is %f',D_a/(boltzmann_const*temperature))
sprintf('Actin forces will have a value %d',P_max)
sprintf('Pinning forces given by %d',pinning_strength)
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

% Variation amplitude (±5% variation)
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


%% ================= TIME EVOLUTION SETUP =================

time_points = (0:num_time_steps-1) * time_step;

% Initialize arrays
height_current = 1e-10 * randn(simulation_grid_size);

height_time_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);

height_time_series(:, :, 1) = height_current;
fprintf('Starting simulation with %d steps...\n', num_time_steps);

% phi initialization (example)
protein_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
% lipid_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
% curvature_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
% kappa_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);
% gamma_map_series = zeros (simulation_grid_size, simulation_grid_size, num_time_steps);
% cytoskeleton_map_series = zeros(simulation_grid_size, simulation_grid_size, num_time_steps);



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
qcut = 0.5*(2*pi/(2*grid_spacing));
filter_q = exp(-(qphys/qcut).^8);   % soft rolloff

figure();
imagesc(qphys);
sprintf('The cutoff is %f',qcut);
figure();
imagesc(filter_q);


% constants for diagnostics or alternative schemes
mobility_det = 1/(4*viscosity * qcut);  % characteristic (not used in modal update)
% End precompute


%% =================  MAIN LOOP  =================
fig1 = figure('Position', [100, 100, 600, 600], 'Color', 'w');
fig2 = figure('Position', [100, 100, 800, 800], 'Color', 'w');


% -----------------------------------------------------------------------
if ~exist('cs_presence','var')
    cs_init_prob = 0.65; % initial fraction of sites with cytoskeleton
    protrusion_init_prob = 0.35; 
    intrusion_init_prob = 0.1;
    cs_presence = double(rand(simulation_grid_size) < cs_init_prob);
    protrusion_map = double(rand(simulation_grid_size) < protrusion_init_prob);
    intrusion_map = double(rand(simulation_grid_size) < intrusion_init_prob);

end
eta_cs = 0.5;
channel_map = zeros(N,N);
total_area_change = 0;
I_grad_all = [];

xi_active = zeros(1, num_time_steps);
% --- 1. Initialize Active Field ---
    Xi_active_q = zeros(N,N); 
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
    
    
    
    psi_binary = double(psi > 0);

%--------------------------------- Activity parameters--------------------------------------


%--------------------------------- Activity parameters--------------------------------------
    
%   1 --------------------------- CYTOSKELETON TURNOVER (stochastic) --------------------------------
    % Use s_activity and tau_cs to set per-step attach/detach probabilities.
    % If currently attached: detach with p_detach = dt / tau_cs
    % If currently detached: attach with p_attach = s_activity * dt / tau_cs
    
    p_detach = ((time_step*step))./ tau_cs;
    p_attach = ((s_activity)*(time_step*step)) ./ tau_cs;
    p_detach = min(p_detach, 1);
    p_attach = min(p_attach, 1);
    
    randmat = rand(simulation_grid_size);
    detach_idx = (cs_presence == 1) & (randmat < p_detach);
    attach_idx = (cs_presence == 0) & (randmat < p_attach);
    cs_presence(detach_idx) = 0;
    cs_presence(attach_idx) = 1;
    pinning_potential = zeros(simulation_grid_size, simulation_grid_size);
    pinned_sites = find(cs_presence>0);
    pinning_potential(pinned_sites) = pinning_strength;

    randmat_protrusion = rand(simulation_grid_size);
    protrude_idx = (cs_presence == 1) & (randmat_protrusion < protrusion_map);
    protrusion_map(protrude_idx) = 1;

    randmat_intrusion = rand(simulation_grid_size);
    intrude_idx = (cs_presence == 1) & (randmat_intrusion < intrusion_map);
    intrusion_map(intrude_idx) = 1;




    % ---- GATE: mobility and reaction rates modulated by cytoskeleton presence ----
    % create gating factor in [0,1] from cs_presence: either binary or smoothed
    cs_gate = cs_presence; % binary gating; replace by smoothed cs_prob_map if you created it

    % ensure mobility is higher where no cytoskeleton is present
    mobility_phi_loc = mobility_phi_max .* (~cs_gate);

    % ---- PROTEIN DYNAMICS AS DENSITY (conserved CH + local production/decay) ---

    % 1) CH-conserved flux term (spectral implementation with spatially uniform mobility assumes constant mobility).
    %    For spatially-varying mobility we do the conservative form in real space:
    %     phi_t_CH = -div( M(x) * grad(mu_phi) )
    %    Implement approximately by computing J = - M .* grad(mu_phi) in real space, then take divergence.
    
    
    coeff1 = 1 ./ (1 - phi) - A_hat .* phi;
    coeff2 = 1 ./ (1 - phi).^2 - A_hat;
    coeff3 = (A_hat ./ (2*S_hat));
    term1 = Laplacian_phi .* coeff1;
    term3 = -phi .* (coeff3 .* Biharmonic_phi);
    term2 = grad_sq_phi .* coeff2;
    
    % Chemical Potential for proteins
    mu_phi = term1 + term3 + term2 - kappa_map .* (Laplacian_h - C_p .* phi);
    mu_phi_q = fft2(mu_phi);
    grad_mu_x = real(ifft2(1i * Qx_p .* mu_phi_q));
    grad_mu_y = real(ifft2(1i * Qy_p .* mu_phi_q));
    % fluxes
    Jx = - mobility_phi_loc .* grad_mu_x;
    Jy = - mobility_phi_loc .* grad_mu_y;
    % divergence of flux (use spectral derivative)
    Jx_q = fft2(Jx);
    Jy_q = fft2(Jy);
    divJ = real(ifft2(1i * Qx_p .* Jx_q + 1i * Qy_p .* Jy_q)); % note: real(ifft2(...)) gives real-space divergence
    % conservative CH update contribution 
    phi_CH_update = time_step .* divJ; 
    
    % 2) Local production/decay: proteins are produced preferentially on rafts AND where cs is present
    % scale k_on by psi_binary (raft) and cs_gate (cytoskeleton)
    k_on_local = k_on_base .* (1 + J .* psi_binary);  % J increases attachment in rafts
    production_term = k_on_local .* cs_gate .* (max_local_density - phi);  % limited by max density
    decay_term = - k_off .* phi .* (1 - cs_gate); % decay faster where no cytoskeleton 
    
    % 3) Optionally add small non-conserved local exchange with cytosol (helps numerics)
    phi_new = phi + phi_CH_update + time_step .* (production_term + decay_term);
    
    % 4) stochastic noise for aggregation only where cs is present (gate the noise)
    % Use your Hermitian noise Xi but multiply by cs_gate in real-space to disable noise where no cs
    % build a real-space mask:
    cs_mask_real = cs_gate;
    % for a spectral noise consistent with earlier stoch_q, transform cs_mask and multiply stoch_q_q by it:
    cs_mask_q = fft2(cs_mask_real);
    % existing stoch part (if you use the prior spectral OU noise), gate it in q-space:
    % (recompute stoch_q as before, then gate:)
    % stoch_q = sqrt( Vq .* (1 - alpha_q.^2) ) .* Xi .* filter_q;
    % stoch_q = stoch_q .* cs_mask_q;  % (this gates modes where cs absent)
    %
    % If you prefer simpler: add small real-space noise only where cs present:
    noise_real = (randn(simulation_grid_size) .* (cs_mask_real)) * 1e-3; % tune amplitude
    phi_new = phi_new + noise_real;
    
    % 5) enforce non-negativity and a reasonable upper bound
    phi_new(phi_new < 0) = 0;
    phi_new(phi_new > max_local_density) = max_local_density;
    
    % ---- LIPID MAP: make binary raft/non-raft explicit (you already did psi_binary) ----
    % If you want psi_new to remain binary for next steps, reassign:
    psi_new_binary = psi_binary;  % keep a binary map for the coupling below
    % But if you want psi to still evolve continuously (and then binarize for readout), keep psi_new as is.
    
    % ---- Mutual preference (optionally bias chemical potentials / kinetics) ----------
    % To bias lipids to become raft where phi large and vice versa, add coupling:
    % Example: increase k_on_local where psi_binary==1 and phi is high (already used J above).
    % Also adjust mu_psi to include -J*(phi - psi) (you had mu_psi = alpha*(psi^3-psi) - lambda*Lap_psi - J*(psi-phi);)
    % That's consistent: larger phi favors psi -> raft formation energetically
    
    %Spatially varying mobility for lipids: higher near cytoskeleton
    mobility_psi_local = mobility_psi .* (1 + eta_cs .* (~cs_presence)); 
    
    % Chemical potential (same form but using latest phi_new)
    mu_psi = alpha .* (psi.^3 - psi) + lambda .* Laplacian_psi - J * (psi - phi_new);
    
    % Compute gradient of mu_psi 
    mu_psi_q = fft2(mu_psi);
    grad_mu_psi_x = real(ifft2(1i * Qx_p .* mu_psi_q));
    grad_mu_psi_y = real(ifft2(1i * Qy_p .* mu_psi_q));
    
    % Compute lipid flux
    Jx_psi = -mobility_psi_local .* grad_mu_psi_x;
    Jy_psi = -mobility_psi_local .* grad_mu_psi_y;
    
    % Divergence of flux
    Jx_psi_q = fft2(Jx_psi);
    Jy_psi_q = fft2(Jy_psi);
    divJ_psi_q = real(ifft2(1i * Qx_p .* Jx_psi_q + 1i * Qy_p .* Jy_psi_q));
    
    % Time update for ψ
    psi_new = psi + time_step .* divJ_psi_q;
    
    % Optional weak cytoskeletal bias toward raft state (stabilize ψ≈1 on cs)
    
    psi_new = psi_new + J_l * cs_presence .*  psi_new;
    
    
    % Enforce ψ bounds to stay physical (0–1 if binary)
    psi_new = min(max(psi_new, 0), 1);
    % Binary raft mask (for coupling with proteins)
    psi_binary = double(psi_new > 0.2);
 
    % ---- STORE updates (replace previous phi/psi updates) -----------------------
    protein_map_series(:,:,step) = phi_new;
    lipid_map_series(:,:,step) = psi_new;
    cytoskeleton_map_series(:,:,step) = cs_presence;
    
    % update kappa_map & C0_map as before (using phi_new and psi_new or psi_binary depending on model)
    
    % --- Candidate kappa (local stiffening/softening) ---
    kappa_candidate = base_kappa .* (1 + beta_phi .* phi_new + beta_psi .* psi_new_binary);
    % clamp to physical bounds
    kappa_candidate(kappa_candidate < base_kappa*0.1) = base_kappa*0.1;
    kappa_candidate(kappa_candidate > base_kappa*10)  = base_kappa*10;
    
    % --- Candidate tension: allow protein to decrease local tension (gamma_phi < 0) ---
    gamma_candidate = base_gamma .* (1 + gamma_phi .* phi_new + gamma_psi .* psi_new_binary);
    
    % --- add curvature-alignment term to further lower tension where curvature matches protein preference ---
    chi = 0.1;  % tune strength
    gamma_candidate = gamma_candidate - chi .* (Laplacian_h - C_p .* phi_new).^2;
    
    % --- enforce global mean tension (area constraint) ---
    gamma_mean_candidate = mean(gamma_candidate(:));
    gamma_map = gamma_candidate .* (base_gamma ./ gamma_mean_candidate);
    
    % --- finally assign kappa_map and gamma_map used in bending/tension forces ---
    kappa_map = kappa_candidate;
    
  
    kappa_map_series(:,:,step) = kappa_map;
    gamma_map_series(:,:,step) = gamma_map;
    C_l = C_ld * ones(size(psi));
    C_l(psi_new_binary > 0) = C_lo;      
    C0_map = C_p.* (phi_new./max_local_density) + C_l;
    
    
    
    
    %   2 --------------------------- VESICULATION EVENTS (stochastic) --------------------------------
    
    % --- normalized tension field (γ-map) ---
    gamma_norm = (gamma_map - min(gamma_map(:))) / (max(gamma_map(:)) - min(gamma_map(:)) + eps);
    
    % --- bias factors ---
    bias_exo  = gamma_norm.^2;        % stronger bias toward high γ
    bias_endo = (1 - gamma_norm).^2;  % stronger bias toward low γ

    
    % --- cytoskeleton mask (0 = free membrane, 1 = anchored) ---
    mask = (cs_presence == 1);             % only allow events on cytoskelton attached membrane
    
    % --- local event probabilities ---
    p_exo_map  = rate_exo_per_m2 * (grid_spacing^2) * time_step * bias_exo .* mask;
    p_endo_map = rate_endo_per_m2 * (grid_spacing^2) * time_step * bias_endo .* mask;
    
    % r0_px = 10; % radius of effect in pixels
    % amp_h_m = 20e-9;     % depth of pit (m) for endocytosis (40 nm)
    % amp_phi = 0.5;       % additional protein density (dimensionless)
    % amp_psi = 0.3;       % additional raft fraction
    % dkappa_amp = 5e-20;  % additional bending rigidity (J)
    % dgamma_amp = 0.0;    % local tension change (N/m)
     
    % --- stochastic selection of protrusion and intrusion from the generated areas of protrusion---
    randmat_protrusion = rand(N, N);
    protrusion_sites  = find(randmat_protrusion<protrusion_map);
    randmat_intrusion = rand(N, N);
    intrusion_sites = find(randmat_intrusion<intrusion_map);
    event_map  = zeros(N,N);
    event_map(protrusion_sites) = 1;
    dx = grid_spacing;
    dy = dx;


    
    % % apply endocytosis events
    % for idx = endo_sites'
    %     [iy, ix] = ind2sub([N,N], idx);
    %     [dh, dphi, dpsi, dkappa, dgamma] = vesicle_kernel(N,N,iy,ix,r0_px, amp_h_m, amp_phi, amp_psi, dkappa_amp, dgamma_amp);
    %     % ramp over N_ramp steps (optional)
    %     h = h + dh;
    %     phi = max(0, phi - dphi); % remove (or concentrate) proteins appropriately
    %     psi = max(0, psi - dpsi);
    %     kappa_map = kappa_map + dkappa;
    %     gamma_map = gamma_map + dgamma;
    %     total_area_change = total_area_change - sum(dh(:))*(dx*dy); % area bookkeeping
    % end
    % 
    % % apply exocytosis events
    % for idx = exo_sites'
    %     [iy, ix] = ind2sub([N,N], idx);
    %     [dh, dphi, dpsi, dkappa, dgamma] = vesicle_kernel(N,N,iy,ix,r0_px, amp_h_m, amp_phi, amp_psi, dkappa_amp, dgamma_amp);
    %     h = h - dh;
    %     phi = max(0, phi - dphi); 
    %     psi = max(0, psi - dpsi);
    %     kappa_map = kappa_map + dkappa;
    %     gamma_map = gamma_map - dgamma;
    %     total_area_change = total_area_change + sum(dh(:))*(dx*dy); % area bookkeeping
    % end
    %
    % sigma_area = K_A *(total_area_change/(membrane_area*1e6));
    % sigma_mean = base_gamma + sigma_area;
    % gamma_map = gamma_map + 0.1 * (sigma_mean - mean(gamma_map(:)));
    
    %Generation of force fields for the membrane
    grad_sq_map = h_x.^2 + h_y.^2;
    I_grad = sum(grad_sq_map(:)) * (grid_spacing^2);
    Sigma_area = (K_A / (2 * membrane_area)) * I_grad;
    I_grad_all = [I_grad_all,I_grad];
    
    % Correct bending force with spatially varying kappa:
    field = kappa_map .* (Laplacian_h - C0_map);
    F_bend_full = -real(ifft2(-Q2 .* fft2(field)));
    
    % Correct tension force with spatially varying gamma:
    
    gx = gamma_map .* h_x;
    gy = gamma_map .* h_y;
    F_tens_full = real(ifft2( 1i*(Qx .* fft2(gx) + Qy .* fft2(gy)) ));  % = ∇·(gamma ∇h)
    
    % Pinning / external (local)
    F_pin = -pinning_potential*height_current;
    
    % Actin protrusion forces
    G_act  = exp(-step*time_step/tau_withdrawal);
    F_protrusion = (P_max*dx*dy)*protrusion_map.*G_act;
    
    
    % ------------------- Compose FULL real-space force ------------------------------
    F_full_real = F_bend_full + F_tens_full;
    
    
    
    % ------------------- Build diagonal Kq using base constants ----------------
    kappa_mean = mean(kappa_map(:));
    gamma_mean = mean(gamma_map(:));
    
    % The spectral modulus is given from the diagonal assumption
    Kq = kappa_mean .* Q4 + gamma_mean .* Q2;
    Kq(1,1) = Inf; % avoid division by zero for q=0
    
    
    % ------------------- Extract residual (non-diagonal) forces ---------------------
    % residual = full force - (linear operator applied in real space)
    % linear_op_real = - base_kappa * ∇^4 h + base_gamma * ∇^2 h + Sigma_area * ∇^2 h + F_press_linear
    % Compute linear_op_real in real space (use same spectral derivatives)
    linear_bend_real = - base_kappa .* Biharmonic_h ;
    linear_tens_real = base_gamma .* Laplacian_h;
    linear_op_real = linear_bend_real + linear_tens_real ;
    
    % Residual forcing (what goes into F_nl_real)
    F_nl_real_residual = F_full_real - linear_op_real;
    % zero out tiny numerical noise
    F_nl_real_residual(abs(F_nl_real_residual) < 1e-20) = 0;
    
    % ------------------- Transform residual into q-space (dealias / mask optionally) -
    Fq_nl = fft2(F_nl_real_residual);
    Fq_nl = Fq_nl .* filter_q;   % if you use a de-alias mask (recommended)
    
    
    
    % ------------------- OU modal parameters ----------------
    tau_q = 1 ./ (Lambda_q .* Kq);
    % --- mobility and Kq already computed earlier ---
    % tau_q = 1 ./ (Lambda_q .* Kq);   % you already have this
    % Compute alpha using time_step (NOT step)
    alpha_q = exp(- time_step ./ tau_q);
    alpha_q(isnan(alpha_q)) = 0;
    
    
    
    
    % 9) deterministic modal forcing from residual
    det_nl_q = (1 - alpha_q) .* (Fq_nl ./ Kq);
    det_nl_q(isnan(det_nl_q) | isinf(det_nl_q)) = 0;
    
    % Generate Hermitian Gaussian noise
    Xi = (randn(simulation_grid_size, simulation_grid_size) + 1i*randn(simulation_grid_size, simulation_grid_size));
    Xi = makeHermitian(Xi);   % ensure conjugate symmetry and Nyquist/DC real
    Xi(1,1) = real(Xi(1,1));
    
    Vq = mobility_det* (boltzmann_const * temperature) ./ (Kq * (grid_spacing^4));
    % clean up numerics
    Vq(isinf(Vq) | isnan(Vq)) = 0;
    % stochastic increment (per-mode) explicit time evolution
    stoch_q = sqrt( Vq .* (1 - alpha_q.^2)) .* Xi.*filter_q;
    
    % 2. Evolve the Active Force Field (Ornstein-Uhlenbeck)
    % We evolve independent noise for each mode (spatially uncorrelated active noise)
    
    % Update prefactor (Variance D_a, Correlation tau_a)
    %      mobility_det* (boltzmann_const * temperature) ./ (Kq * (grid_spacing^4))
    ou_prefactor = sqrt(2 * D_a^2 * tau_a)./ (Kq * (grid_spacing^4));
    
    % Generate fresh white noise for every mode
    eta_q = (randn(size(Hq)) + 1i*randn(size(Hq)));
    eta_q = makeHermitian(eta_q); % Ensure real output in real-space
    % (Note: Scale eta_q by 1/dx if needed for unit consistency, usually implied in D_a)

    % Euler-Maruyama update for the Force Field
    d_Xi = -(Xi_active_q ./ tau_a) * time_step + ou_prefactor .* eta_q;
    Xi_active_q = Xi_active_q + d_Xi + F_pin + F_protrusion;

    % 3. Apply as a Force to the Membrane
    % This active force acts like the "Nonlinear/Explicit" part.

    Active_Kick = Lambda_q .* Xi_active_q * time_step;

    % If using the ETD1 scheme (your previous code), it's treated like the deterministic force:
    % The term (1-alpha)/Kq is effectively (Mobility * dt) for high frequencies.
    % So rigorous integration into ETD1:
    Active_Term = ((1 - alpha_q)./Kq) .* (Xi_active_q); 
    % Note: Xi is treated as Force here. (Force / Stiffness = Displacement)
    
    % 11) final modal update 
    % Update OU Process
    Hq_new = alpha_q .* Hq + det_nl_q + stoch_q + Active_Term;
    
    % 12) optional spectral filter and back to real
    Hq_new = Hq_new .* filter_q;
    height_current = real(ifft2(Hq_new));
    
    % 13) remove mean (control q=0 / volume)
    height_current = height_current - mean(height_current(:));
    
    % 14) store & visualize
    height_time_series(:,:,step) = height_current;
    
    
    
    
    % Update visualization less frequently
    if mod(step,10000) == 1
        update_visualization(fig1,fig2, X_nm, Y_nm,height_current, time_points(step),lipid_map_series,kappa_map_series,protein_map_series,step,cs_presence,boltzmann_const,temperature);
        drawnow;
    end
end

%% Computing time averaged PSD of entire simulation stack 



fprintf('Computing normalized global temporal PSD across all q modes...\n');

dx_nm = grid_spacing*1e9;     % nm
dy_nm = dx_nm;
[Ny, Nx, Nt] = size(height_time_series);
dt_s = time_step;

% Convert height to nm
h_nm_series = height_time_series * 1e9;

% Preallocate global height fluctuation signal (units: nm²)
global_time_series = zeros(Nt, 1);

for t = 1:Nt
    h_frame = h_nm_series(:,:,t);

    % remove mean height (DC offset)
    h_frame = h_frame - mean(h_frame(:));

    % compute spatially averaged height-squared (nm²)
    global_time_series(t) = mean(h_frame(:).^2);
end

% --------------- Temporal PSD (nm²/Hz) ------------------

nfft_global = 2^nextpow2(Nt);
window = hann(Nt)';
sig_win = global_time_series' .* window;

spec_global = fft(sig_win, nfft_global);

freq_Hz_global = (0:nfft_global/2-1) / (nfft_global * dt_s);

% PSD in nm² / Hz
PSD_global_norm = (2 * dt_s / (Nt * sum(window.^2))) ...
                   * abs(spec_global(1:nfft_global/2)).^2;


% 9. Theoretical PSD(f) using fluctuation spectrum integral
% -------------------------------------------------------------------------

% Physical constants
kB = 1.380649e-23;       % J/K
T_kelvin = 300;          % temperature
Activity = 0;    % membrane area in m^2

% Use mean bending rigidity and tension from maps
kappa_mean = mean(kappa_map(:));           
sigma_mean = mean(gamma_map(:));         
gamma_local = 0;                 
mean_confinement = pinning_strength;
% Effective viscosity (your η_eff)
eta_eff = viscosity;     

% Membrane shear modulus μ
mu = 0;               

% q-range from grid spacing
Lx = simulation_grid_size * grid_spacing; 
Ly = simulation_grid_size * grid_spacing;

qmin = 2*pi / max(Lx, Ly);
qmax = qcut;       

Nq = 1000;
q = linspace(qmin, qmax, Nq);

% Preallocate theoretical PSD
PSD_theoretical = zeros(size(freq_Hz_global));

% Main loop over frequencies
for i = 1:length(freq_Hz_global)

    f = freq_Hz_global(i);
    omega = 2*pi*f;
    
    
    % Frequency term
    denom_freq = (4 * eta_eff * (2*pi*f)).^2;
    
    % Present theory in  IRM paper
    denom_mech_1 = ( ...
        kappa_mean .* q.^3 + ...
        sigma_mean .* q + mean_confinement./q).^2;
    
    % Considering only thermal fluctuations
    denom_mech_2 = ( ...
        kappa_mean .* q.^3 + ...
        sigma_mean .* q.^1 + mean_confinement./q).^2;
        
  
    % ---------- active contribution ----------
    % choose spatial spectrum of active forcing S_Fq(omega) = 2 D_a(q) / (1 + (omega tau_a)^2)
    % We'll set D_a(q) = Da0 * exp(- (q/qc)^2) as example; choose qc if needed
    qc = qcut;  % 1/m example length scale (tune)
    Da_q = D_a .* exp( - (q./qc).^2 );  
    
    % integrand_1  is for IRM paper
    integrand_1 = 1 ./ (denom_freq + denom_mech_1);
    %integrand_2 is for pure thermal fluctuations
    integrand_2 =  1./((omega)^2 + denom_mech_2);
    
    %integrand_3 is for frequency dependent activity along with thermal
    %fluctuations
    S_F_q_omega = 2*Da_q./ (1 + (omega * tau_a).^2);
    integrand_3  = 1 ./ ((omega)^2 + denom_mech_2) .* S_F_q_omega; 
    
    % active contribution to mode amplitude PSD: |G|^2 * S_F  where |G|^2 = 1/(omega^2 + omega_q^2)
    % but in our mechanical mapping |G|^2 corresponds to 1/denom (we already have denom)
    % scale prefactor for area similar to thermal prefactor: active prefactor should include geometry
    
    PSD_theoretical_irm(i) = (4*eta_eff*D_a/pi) * trapz(q, integrand_1);
    
    PSD_theoretical(i) = ((4 * eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_2) + ((4 * eta_eff)/ pi) * trapz(q,integrand_3);
    
    PSD_theoretical_passive(i) = ((4 * eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_2);
    
end

% convert from m^2/Hz → nm^2/Hz
PSD_theoretical = PSD_theoretical * 1e18;
PSD_theoretical_passive = PSD_theoretical_passive*1e18;
PSD_theoretical_irm = PSD_theoretical_irm * 1e18;
% -------------------------------------------------------------------------
% Plot measured vs theoretical PSD
% -------------------------------------------------------------------------


figure('Color','w'); 
loglog(freq_Hz_global, PSD_global_norm, 'g--', 'LineWidth', 1); hold on;
loglog(freq_Hz_global, PSD_theoretical, 'r--', 'LineWidth', 2);
loglog(freq_Hz_global, PSD_theoretical_irm, 'b--', 'LineWidth', 2);
loglog(freq_Hz_global, PSD_theoretical_passive, 'm--', 'LineWidth', 2);



xlabel('Frequency [Hz]');
ylabel('Power [nm^2/Hz]');
legend('Measured PSD', 'New Theoretical PSD','IRM theoretical PSD','Passive PSD','Location','SouthWest');
title('Measured vs Theoretical PSD of Membrane Fluctuations');
grid on;

%%
% cd(dir);
% % Setup for video
% video_filename = 'membrane_dynamics_movie(pinned).mp4';
% video1 = VideoWriter(video_filename, 'MPEG-4');
% video1.FrameRate = 100;
% open(video1);
% figure();  % Create and reuse figure window
% h_range = max(abs(height_time_series(:)));
% 
% cmin = -h_range*1e9;
% cmax=  h_range*1e9;
% for t = 1:num_time_steps/2
%     clf;
%     h = height_time_series(:,:,t)*1e9;
%     surf(X_nm, Y_nm, h, 'EdgeColor','none');
%     hold on;
% 
%     % Use actual time value in microseconds
%     title(sprintf('Membrane Height (t = %.2f ms)', time_points(t)*1e3));
%     xlabel('x/L'); ylabel('y'); zlabel('h (nm)');
%     view(-30, 60);
%     colormap(turbo);
%     hcb = colorbar('Location','eastoutside');
%     hcb.Label.String = 'Height (nm)';
%     caxis([cmin cmax]);
%     zlim([cmin cmax]);
%     light('Position', [1 1 1], 'Style', 'infinite');
% 
%     drawnow;
%     frame = getframe(gcf);
% 
%     writeVideo(video1, frame);
% end
% 
% close(video1);
% cd(old_dir);

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

% for step = 1:50:num_time_steps
%     
%         height_current = height_time_series(:,:,step);
%        
%         phi_norm = mat2gray(protein_map_series(:,:,step));
%         phi_image = im2uint8(phi_norm);
%         filename = sprintf('Protein_Density_Map_%04d.png',step);
%         imwrite(phi_image,fullfile(output_dir_protein_map,filename));
%         
%         
%         a = lipid_map_series(:,:,step);
%         binary_lipid_map_series(:,:,step) = a > 0.01;
%         psi_norm = mat2gray(lipid_map_series(:,:,step));
%         phi_image = im2uint8(psi_norm);
%         filename = sprintf('Lipid_Density_Map_%04d.png',step);
%         imwrite(phi_image,fullfile(output_dir_lipid_map,filename));
%         
%         height_norm = mat2gray(height_time_series(:,:,step)*1e9);
%         height_image = im2uint8(height_norm);
%         filename = sprintf('Height_Map_%04d.png',step);
%         imwrite(height_image,fullfile(output_dir_curvature_map,filename));
%         
%         cyto_norm = mat2gray(cytoskeleton_map_series(:,:,step));
%         cyto_image = im2uint8(cyto_norm);
%         filename = sprintf('Cyto_Map_%04d.png',step);
%         imwrite(cyto_image,fullfile(output_dir_cyto_map,filename));
%         
%         bend_rigid_map = mat2gray(kappa_map_series(:,:,step));
%         bend_rigid_map_image = im2uint8(bend_rigid_map);
%         filename = sprintf('BendRidigity_Map_%04d.png',step);
%         imwrite(bend_rigid_map_image,fullfile(output_dir_bending_map,filename));
%         
%         tension_map = mat2gray(gamma_map_series(:,:,step));
%         tension_image = im2uint8(tension_map);
%         filename = sprintf('Tension_Map_%04d.png',step);
%         imwrite(tension_image,fullfile(output_dir_tension_map,filename));
% 
% end

%% Ensuring all the parameters of the succesful simulation is saved 

readme_file = fullfile(dir,'README.txt');
fid = fopen(readme_file, 'w');

% Write header
fprintf(fid, '=============================================\n');
fprintf(fid, '   Simulation Parameters Log File (README)\n');
fprintf(fid, '=============================================\n\n');

% Example metadata
fprintf(fid, 'Date & Time: %s\n', datestr(now));
fprintf(fid, '---------------------------------------------\n');
fprintf(fid, 'Grid size (NxN): %d x %d\n', N, N);
fprintf(fid, 'Grid spacing (m): %.3e\n', grid_spacing);
fprintf(fid, 'Time step (s): %.3e\n', time_step);
fprintf(fid, 'Number of steps: %d\n\n', num_time_steps);

% Write physical parameters
fprintf(fid, '--- Physical Constants ---\n');
fprintf(fid, 'k_B T (J): %.3e\n', boltzmann_const * temperature);
fprintf(fid, 'Base bending modulus (kappa): %.3e\n', base_kappa);
fprintf(fid, 'Base tension (gamma): %.3e\n', base_gamma);

% Write protein/lipid parameters
fprintf(fid, '--- Protein / Lipid Coupling ---\n');
fprintf(fid, 'Mobility_phi: %.3e\n', mobility_phi);
fprintf(fid, 'Mobility_psi: %.3e\n', mobility_psi);
fprintf(fid, 'Coupling J: %.3e\n', J);
fprintf(fid, 'C_p (protein curvature): %.3e\n', C_p);
fprintf(fid, 'C_l (lipid curvature): %.3e\n', C_l);
fprintf(fid, 'beta_phi (protein stiffening): %.3e\n', beta_phi);
fprintf(fid, 'beta_psi (lipid stiffening): %.3e\n\n', beta_psi);

% Cytoskeleton / activity section
fprintf(fid, '--- Cytoskeleton Turnover ---\n');
fprintf(fid, 's_activity: %.3e\n', s_activity);
fprintf(fid, 'tau_cs (s): %.3e\n', tau_cs);
fprintf(fid, 'Initial CS coverage: %.2f\n', mean(cs_presence(:)));

% Derived fields
fprintf(fid, '--- Derived Quantities ---\n');
fprintf(fid, 'Lambda_q range: [%g, %g]\n', min(Lambda_q(:)), max(Lambda_q(:)));
fprintf(fid, 'Filter_q applied: %d\n', any(filter_q(:)~=1));
fprintf(fid, '\n');

% Close file
fclose(fid);

fprintf('✅ README.txt created at: %s\n', readme_file);

%% ================= VISUALIZATION FUNCTION =================
function [dh, dphi, dpsi, dkappa, dgamma] = vesicle_kernel(Ny, Nx, iy, ix, r0, amp_h, amp_phi, amp_psi, dkappa_amp, dgamma_amp)
    % create Gaussian kernel on grid
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    R2 = (X-ix).^2 + (Y-iy).^2;
    kernel = exp(-R2/(2*(r0^2)));
    kernel = kernel / sum(kernel(:));   % normalize to unit integrated weight

    % dh: change in height (nm) or meters depending on units. For budding (endo) we lower center -> negative dh
    dh = amp_h * kernel;   % negative -> invagination; positive -> bulge for exocytosis
    dphi = amp_phi * kernel; % concentrate proteins into pit
    dpsi = amp_psi * kernel; % concentrate raft lipids if applicable
    dkappa = dkappa_amp * kernel;  % local stiffening
    dgamma = dgamma_amp * kernel;  % local change in membrane tension (could be negative for exo)
end

function update_visualization(fig1,fig2,X_nm, Y_nm, h_total, current_time,lipid_map_series,kappa_map_series,protein_map_series,step,cs_presence,boltzmann_const,temperature)
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
    subplot(2,2,1);
    contourf(X_nm, Y_nm, protein_map_series(:,:,step), 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Protein map');


    subplot(2,2,2);
    contourf(X_nm, Y_nm, kappa_map_series(:,:,step)/(boltzmann_const*temperature), 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Bending rigidity map');
    
    a = lipid_map_series(:,:,step);
    subplot(2,2,3);
    contourf(X_nm, Y_nm, a, 30, 'LineColor','none'); % 30 contour levels
    axis image; colorbar;
    colormap(turbo);
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Lipid Phase Map');
    drawnow;
    
    subplot(2,2,4);
    imshow(cs_presence,[]);
    axis image;
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Cytoskeleton Map');
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

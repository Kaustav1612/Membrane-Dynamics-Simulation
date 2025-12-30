% Membrane Dynamics Simulation with Physical Mode Evolution and Advanced Visualization
clear;close all; clc;

%% ================= PHYSICAL PARAMETERS =================

simulation_grid_size = 60;       % Reduced for faster testing
simulated_area =5e-12;        % 5μm x 5μm membrane (15 μm²)
membrane_area = 5e-12;            % Membrane area should match simulated area
viscosity = 10;                 % Water viscosity (Pa·s)
temperature = 310;                % Temperature (K)
boltzmann_const = 1.38e-23;       % Boltzmann constant (J/K)
bending_rigidity = 12 * boltzmann_const * temperature; % ~20 kT
time_step = 10e-6;                 % MUCH smaller time step (20 ms)
total_simulation_time = 5;     % 50 ms total timefTpo
surface_tension = 40e-6;           % MUCH lower tension (40 pN/um)
old_dir = pwd;
dir = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_absent\heterogenous\03(low_tension)';
output_dir_curvature_map = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_absent\heterogenous\03(low_tension)\results_curvature_map';
pinning_strength = 80e-6;          % Moderate pinning (N/m²)
low_kappa = 10;
high_kappa = 15;
low_tension = 50e-6;
high_tension  = 100e-6;
K_A=0;

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
kappa_variation = 0.1 * base_kappa;
gamma_variation = 0.1 * base_gamma;


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


% Soft harmonic confinement at edges

[X_grid, Y_grid] = meshgrid(1:simulation_grid_size, 1:simulation_grid_size);


%% ================= TIME EVOLUTION SETUP =================
num_time_steps = round(total_simulation_time / time_step);
time_points = (0:num_time_steps-1) * time_step;


% Define the amplitude of the sine wave (e.g., 5 nm)
A =0; 

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
mobility_det = 1/(4*viscosity * qcut);  % characteristic (not used in modal update)
% End precompute

%% =================  MAIN LOOP  =================
fig1 = figure('Position', [100, 100, 600, 800], 'Color', 'w');
I_grad_all = [];
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
    tmp1 = kappa_map .* Laplacian_h;             
    %+  real(ifft2( -Q2 .* fft2(tmp2) )); 
    F_bend_full = - real(ifft2( -Q2 .* fft2(tmp1) )) ;

    % Correct tension force with spatially varying gamma:
    gx = gamma_map .* h_x; gy = gamma_map .* h_y;
    F_tens_full = real(ifft2( 1i*(Qx .* fft2(gx) + Qy .* fft2(gy)) ));

    % Otherwise compute global Sigma_area and treat as linear extra tension (goes into Kq)
    grad_sq_map = h_x.^2 + h_y.^2;
    I_grad = sum(grad_sq_map(:)) * (grid_spacing^2);
    Sigma_area = (K_A / (2 * membrane_area)) * I_grad;
    I_grad_all = [I_grad_all,I_grad];



    % ------------------- Compose FULL real-space force ------------------------------
    F_full_real = F_bend_full + F_tens_full;
    % Note: Sigma_area * Lap_h is included here only if you want it also to appear as part of F_nl.
    % But we'll include Sigma_area in Kq (preferred), so subtract its linear part below.

    % ------------------- Build linear diagonal Kq using base constants ----------------
    Kq = kappa_map .* Q4 + (gamma_map).* Q2;
    Kq(1,1) = Inf; 

    % ------------------- Extract residual (non-diagonal) forcing ---------------------
    % residual = full force - (linear operator applied in real space)
    % linear_op_real = - base_kappa * ∇^4 h + base_gamma * ∇^2 h + Sigma_area * ∇^2 h + F_press_linear
    % Compute linear_op_real in real space (use same spectral derivatives)
    linear_bend_real = - base_kappa .* Biharmonic_h ; 
    linear_tens_real = base_gamma .* Laplacian_h;
    linear_op_real = linear_bend_real + linear_tens_real;

    % Residual forcing (what goes into F_nl_real)
    F_nl_real_residual = F_full_real - linear_op_real;

    % Optional: zero out tiny numerical noise
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
    Vq = (boltzmann_const * temperature*mobility_det) ./ (Kq * (grid_spacing^4));
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
    if rem(step,10000)==1
        update_visualization(fig1, X_nm, Y_nm,L_nm,height_time_series(:,:,step), time_points(step),step);
    end   
end



std_dev_time = std(height_time_series,0,3);
figure;
imshow(std_dev_time*1e9,[min(std_dev_time(:)*1e9) max(std_dev_time(:)*1e9)]);
colormap jet;
colorbar;
caxis([min(std_dev_time(:)*1e9) max(std_dev_time(:))*1e9]);
title('Standard Deviation of height(nm)');

for  step = 1 : 50 : num_time_steps
    height_norm = mat2gray(height_time_series(:,:,step));
    height_image = im2uint8(height_norm);
    filename = sprintf('Height_Map_%04d.png',step);
    imwrite(height_image,fullfile(output_dir_curvature_map,filename));
end

%%
% cd(dir);
% % Setup for video
% video_filename = 'membrane_dynamics_movie.mp4';
% video1 = VideoWriter(video_filename, 'MPEG-4');
% video1.FrameRate = 100;
% open(video1);
% figure();  % Create and reuse figure window
% h_range = max(abs(height_time_series(:)));
% 
% cmin = -h_range*1e9 ;
% cmax=  h_range*1e9 ;
% for t = 1:50:num_time_steps
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
Activity = 1e-4;    % activity

% Use mean bending rigidity and tension from maps
kappa_mean = mean(kappa_map(:));           
sigma_mean = mean(gamma_map(:));         
gamma_local = 0;                 
mean_confinement = 0;
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
    
    % Frequency term
    denom_freq = (4 * eta_eff * (2*pi*f)).^2;

    % Mechanical restoring term
    denom_mech = ( ...
        kappa_mean .* q.^3 + ...
        sigma_mean .* q + ...
        (9*kB*T_kelvin*mu*q)/(16*pi*kappa_mean)+...
        (mean_confinement) ./ q +...
        ((K_A .*q)/4)*(mean(I_grad_all))...
    ).^2;

    integrand = 1 ./ (denom_freq + denom_mech);

    PSD_theoretical(i) = (4*eta_eff*kB*T_kelvin/pi) * trapz(q, integrand);
end
% convert from m^2/Hz → nm^2/Hz
PSD_theoretical = PSD_theoretical * 1e18;

% -------------------------------------------------------------------------
% Plot measured vs theoretical PSD
% -------------------------------------------------------------------------

figure('Color','w'); 
loglog(freq_Hz_global, PSD_global_norm, 'b', 'LineWidth', 1.5); hold on;
loglog(freq_Hz_global, PSD_theoretical, 'r--', 'LineWidth', 2);
ylim([1e-20 1e5]);
xlabel('Frequency [Hz]');
ylabel('Power [nm^2/Hz]');
legend('Measured PSD','Theoretical PSD','Location','SouthWest');
title('Measured vs Theoretical PSD of Membrane Fluctuations');
grid on;


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

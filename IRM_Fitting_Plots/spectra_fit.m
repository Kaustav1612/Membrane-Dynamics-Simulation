close all;clc;
I_min    = 912;
k_factor = conv;


%%
% ---- File Import ----
old_dir = pwd;
olddir = cd('F:\Kaustav\IRM\200126\Bleb\002\set');
directory = 'F:\Kaustav\IRM\200126\Bleb\002';
disp('Scanning TIFF files...');
kB_T = 1.38e-23 * 310;
T_kelvin = 310;
kB = 1.38e-23;
grid_spacing = 72e-9;
side_length = FBR(1, 3);
qmin = 2*pi / (side_length * grid_spacing);
qmax = 2*pi / (grid_spacing);

file_list = dir('*.tif'); % Renamed from A to avoid conflict
num_frames = size(file_list, 1);

if num_frames == 0
    error('No TIFF files found.');
end

% Load Stack
image_stack = [];
info = imfinfo(file_list(1).name);
H = info.Height; W = info.Width;
image_stack = zeros(H, W, num_frames, 'uint16');

for k = 1:num_frames
    image_stack(:,:,k) = imread(file_list(k).name);
end
first_frame = image_stack(:,:,1);
hold on;
continue_drawing = true;

% ---- Background Selection ----
while continue_drawing
    imshow(first_frame,[]);  % Display the black and white image
    title('Draw the Background');
    h = drawrectangle('InteractionsAllowed', 'all', 'Color', 'r');
    
    % Extract top-left corner and dimensions
    top_left = h.Position(1:2);
    width = h.Position(3);
    height = h.Position(4);
    
    % Ensure the shape is a square
    side_length = min(width, height);
    square_coords = [top_left;
        top_left + [side_length, 0];
        top_left + [side_length, side_length];
        top_left + [0, side_length];
        top_left];  % Close the square
    
    
    % Confirm the ROI or redraw
    choice = questdlg('Do you want to keep this ROI or redraw?', ...
        'Confirm ROI', ...
        'Keep', 'Redraw', 'Keep');
    if strcmp(choice, 'Redraw')
        delete(h); % Remove the drawn rectangle
        delete(square_plot); % Remove the displayed red square
        continue;  % Restart the loop
    else
        continue_drawing=false;
    end
end
close;

max_row = round(square_coords(3, 2));
min_row = round(square_coords(2, 2));
max_col = round(square_coords(2, 1));
min_col = round(square_coords(1, 1));


bg_img = first_frame(min_row:max_row, min_col:max_col);
simulation_grid_size_row = max_row - min_row;
simulation_grid_size_col = max_col - min_col;

x_coords = (0:simulation_grid_size_col-1) * grid_spacing;
y_coords = (0:simulation_grid_size_row-1) * grid_spacing;
[grid_X, grid_Y] = meshgrid(x_coords, y_coords);
X_nm = grid_X ;
Y_nm = grid_Y ;
bg_height_stack = zeros(simulation_grid_size_row+1,simulation_grid_size_col+1,num_frames);

for k = 1:num_frames
    frame_idx = k;
    
    % ---- B. Intensity → height ----
    intensity_map = double(bg_img);
    height_map = (intensity_map - I_min) / k_factor;
    height_map(height_map < 0) = 0;
    bg_height_stack(:,:,k) = height_map;
end

[Ny, Nx, Nt] = size(bg_height_stack);
dt_s = 50e-3; %(in s)
% 1. Correct Initialization
% We need to know the length of the PSD output first
nfft_global = 2^nextpow2(Nt);
fs = 1/dt_s;
[~, f_dummy] = pwelch(zeros(Nt,1), hann(floor(Nt/2)), [], nfft_global, fs);
num_freq_points = length(f_dummy);

% Preallocate: Rows = grid points, Cols = frequency points
all_PSD = NaN(simulation_grid_size_row+1,simulation_grid_size_col+1, num_freq_points);

fprintf('Starting ensemble PSD averaging of background...\n');

for k = 1:simulation_grid_size_row+1
    for j = 1 : simulation_grid_size_col+1
        h_patch = bg_height_stack(k, j, :);
        h_series = squeeze(h_patch);
        
        % 3. Pre-processing: Remove DC offset
        h_series = abs(h_series - mean(abs(h_series)));
        
        % 4. Compute PSD using pwelch (Units: nm^2/Hz)
        [psd_tmp, freq_Hz_global] = pwelch(h_series, ...
            hann(floor(Nt)), ...
            [], ...
            nfft_global, ...
            fs);
        
        % Store in matrix
        all_PSD(k,j,:) = psd_tmp';
    end
    
end

% 5. Ensemble Average across all spatial points squeeze converts the 1x1xN result into an Nx1 vector
averaged_bg_PSD = squeeze(mean(all_PSD, [1, 2], 'omitnan'));
freq_Hz_global = freq_Hz_global(:);
cd(old_dir);
%%
% Calculation of Tension and Activity from ensemble height fluctutations and excess area and correlation time

map_sigma = zeros(size(FBR,1), 1);
map_activity = zeros(size(FBR,1), 1);
num_rois = size(FBR,1);



% Tension from excess area calculation for FBRs
% Tension from excess area calculation for FBRs
pixel_size = grid_spacing;
kappa_mean = 12*kB_T;
sigma_excess_area = zeros(2, 1);

% --- 0. Setup the qcut Sweep Array ---
num_q_steps = 100;

% Define fixed qmax (bending-dominated regime limit)
qmax_fixed = 3.5e7; % m^-1 - optical diffraction limit/bending regime boundary

% Create qcut sweep from low to high, avoiding qmax
qmin = 2*pi / (side_length * pixel_size); % Minimum q from system size
qcut_sweep = linspace(qmin, qmax_fixed * 0.95, num_q_steps); % Stop at 95% of qmax

fprintf('Starting q_cut sweep across %d ROIs...\n', size(1,1));
sigma_sweep_all = zeros(size(1,1), num_q_steps);

for roi_idx = 1:5:size(FBR,1)
    fprintf('Processing ROI %d of %d...\n', roi_idx, size(FBR,1));
    
    % --- 1. Extract Height Map ---
    x = FBR(roi_idx, 1);
    y = FBR(roi_idx, 2);
    side_length = FBR(roi_idx, 3);
    
    roi_stack = image_stack(y:y+side_length, x:x+side_length, :);
    h_stack_nm = (double(roi_stack) - I_min) / k_factor;
    h_stack_m = h_stack_nm*1e-9; % Convert to meters
    
    % --- 2. Remove mean and smooth ---
    h_stack_smooth = h_stack_m - mean(h_stack_m(:));
    h_stack_smooth = imgaussfilt(h_stack_smooth, 1.0);
    
    % --- 3. Prepare Fourier Space (q-grid) ---
    [Ny, Nx, Nt] = size(h_stack_smooth);
    dqx = 2*pi / (Nx * pixel_size);
    dqy = 2*pi / (Ny * pixel_size);
    
    qx = ( (0:Nx-1) - floor(Nx/2) ) * dqx;
    qy = ( (0:Ny-1) - floor(Ny/2) ) * dqy;
    [QX, QY] = meshgrid(qx, qy);
    Q_mag = sqrt(QX.^2 + QY.^2);
    
    % Compute FFT of the entire smoothed stack
    H_stack_q = zeros(size(h_stack_smooth));
    for t = 1:Nt
        H_stack_q(:,:,t) = fftshift(fft2(h_stack_smooth(:,:,t)));
    end
    
    % --- 4. Sweep qcut and calculate tension ---
    term_pre = (kB_T) / (8 * pi * kappa_mean);
    A_const = kappa_mean * qmax_fixed^2;
    
    for q_idx = 1:num_q_steps
        current_qcut = qcut_sweep(q_idx);
        B_const = kappa_mean * current_qcut^2;
        
        % Create bandpass mask (keep modes between current_qcut and qmax_fixed)
        q_mask = (Q_mag >= current_qcut) & (Q_mag <= qmax_fixed);
        
        % Apply mask and reconstruct filtered images
        filtered_h_stack = zeros(size(h_stack_smooth));
        for t = 1:Nt
            H_q_filtered = H_stack_q(:,:,t) .* q_mask;
            filtered_h_stack(:,:,t) = real(ifft2(ifftshift(H_q_filtered)));
        end
        
        % Calculate dynamic excess area (alpha)
        [dh_dx, dh_dy] = gradient(filtered_h_stack, pixel_size);
        grad_sq_mag = dh_dx.^2 + dh_dy.^2;
        alpha_dynamic = 0.5 * mean(grad_sq_mag(:));
        
        % Solve for tension
        tension_solver = @(sig) (term_pre * log((sig + A_const) ./ (sig + B_const))) - alpha_dynamic;
        
        try
            % Use a wider initial bracket for better convergence
            sigma_sol = fzero(tension_solver, [1e-8, 1e-3]); % Broader range: 1e-8 to 1e-3 N/m
            sigma_sweep_all(roi_idx, q_idx) = sigma_sol;
        catch
            sigma_sweep_all(roi_idx, q_idx) = NaN;
        end
    end
end

disp('Sweep complete. Analyzing results...');

% --- 5. Analyze sweep results ---
% Calculate mean and standard deviation across ROIs
mean_sigma = mean(sigma_sweep_all, 1, 'omitnan');
std_sigma = std(sigma_sweep_all, 0, 1, 'omitnan');

% Define validity mask
valid_mask = (mean_sigma > 1e-8) & (mean_sigma < 1e-3) & ... % Physical tension range
    (qcut_sweep <= qmax_fixed * 0.99) & ... % Below qmax
    (~isnan(mean_sigma));

valid_qcuts = qcut_sweep(valid_mask);
valid_sigmas = mean_sigma(valid_mask);
valid_stds = std_sigma(valid_mask);

if isempty(valid_sigmas)
    error('No valid tension values found. Check your data and parameters.');
end

% --- 6. Plateau detection with improved algorithm ---
% Find the maximum tension in the valid range
[max_tension, max_idx] = max(valid_sigmas);

% Define plateau threshold (adjustable based on data quality)
plateau_threshold = 0.95 * max_tension;

% Find the first index where tension reaches the plateau
% Start from the maximum and work backwards to find the onset
plateau_start_idx = find(valid_sigmas >= plateau_threshold, 1, 'first');

% If plateau detection fails, use a derivative-based method
if isempty(plateau_start_idx)
    % Calculate derivative of tension with respect to q
    dsigma_dq = diff(valid_sigmas) ./ diff(valid_qcuts);
    dsigma_dq_smooth = smoothdata(dsigma_dq, 'movmean', 5);
    
    % Find where derivative becomes small (approaches plateau)
    small_deriv_threshold = 0.1 * max(abs(dsigma_dq_smooth(1:end/3)));
    plateau_start_idx = find(abs(dsigma_dq_smooth) < small_deriv_threshold, 1, 'first');
    
    if isempty(plateau_start_idx)
        % Fallback: use the point before the maximum starts decreasing
        [~, max_idx_local] = max(valid_sigmas);
        plateau_start_idx = max_idx_local;
    end
end

true_tension = valid_sigmas(plateau_start_idx);
best_qcut = valid_qcuts(plateau_start_idx);
tension_error = valid_stds(plateau_start_idx);

% Calculate critical length scale
critical_length_scale_nm = (2 * pi / best_qcut) * 1e9;

% --- 7. Visualize results ---
figure('Position', [100, 100, 800, 600]);

% Plot mean tension with error bars
errorbar(valid_qcuts, valid_sigmas, valid_stds, 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'b');
hold on;

% Mark the detected plateau onset
plot(best_qcut, true_tension, 'ro', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'LineWidth', 2);
xline(best_qcut, 'r--', 'Detected Plateau Onset', ...
    'LabelOrientation', 'horizontal', 'LineWidth', 1.5);

% Add horizontal line at plateau threshold
yline(plateau_threshold, 'g--', 'Plateau Threshold', ...
    'LabelOrientation', 'horizontal', 'LineWidth', 1);

% Format plot
xlabel('q-cut (m^{-1})', 'FontSize', 12);
ylabel('Apparent Tension (N/m)', 'FontSize', 12);
title('q-Cut Sweep for Membrane Tension Estimation', 'FontSize', 14);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([min(valid_qcuts), max(valid_qcuts)]);



hold off;

% --- 8. Output results ---
fprintf('\n========== RESULTS ==========\n');
fprintf('Selected q_cut = %.2e m^-1\n', best_qcut);
fprintf('Critical Thermal Length Scale: %.1f nm\n', critical_length_scale_nm);
fprintf('Estimated True Mechanical Tension: %.2e ± %.2e N/m\n', ...
    true_tension, tension_error);
fprintf('=============================\n\n');

% Optional: Save results
% save('tension_results.mat', 'true_tension', 'tension_error', 'best_qcut', ...
%      'critical_length_scale_nm', 'valid_qcuts', 'valid_sigmas');

qcut = best_qcut;
for roi_idx = 1:size(FBR,1)
    x = FBR(roi_idx, 1);
    y = FBR(roi_idx, 2);
    side_length = FBR(roi_idx, 3);
    roi_stack = image_stack(y:y+side_length, x:x+side_length, :);
    h_stack_nm = (double(roi_stack) - I_min) / k_factor;
    h_stack_m = h_stack_nm * 1e-9; % Convert to meters for physics formulas!
    
    [dh_dx, dh_dy] = gradient(h_stack_m, pixel_size);
    grad_sq_mag = dh_dx.^2 + dh_dy.^2;
    alpha = 0.5 * mean(grad_sq_mag(:));
    
    % 5. Solve for Tension (sigma)
    
    % The equation is non-linear: alpha = C * ln( (sigma + A)/(sigma + B) )
    term_pre = (kB_T) / (8 * pi * kappa_mean);
    A_const = kappa_mean * qmax^2;
    B_const = kappa_mean * qcut^2;
    
    % Define function to find root (Equation - Alpha = 0)
    tension_solver = @(sig) (term_pre * log((sig + A_const) ./ (sig + B_const))) - alpha;
    
    try
        sigma_sol = fzero(tension_solver, [1e-8, 5000e-6]);
        sigma_excess_area(roi_idx) = sigma_sol;
    catch
        sigma_excess_area(roi_idx) = NaN; % Failed to converge (or negative tension)
    end
end

% Calculation of Tension and Activity Map for FBRs

for roi_idx = 1:size(FBR,1)
    
    % --- 2. Calculate the Geometric Factor ---
    kappa_mean = 12*kB_T;
    sigma_mean = sigma_excess_area(roi_idx);
    mean_confinement = 0;
    x = FBR(roi_idx, 1);
    y = FBR(roi_idx, 2);
    
    roi_stack = image_stack(y:y+side_length, x:x+side_length, :);
    h_stack = (double(roi_stack) - I_min) / k_factor;
    h_stack = h_stack .* 1e-9;
    
    % 1. Remove the static cell shape (DC offset for each pixel over time)
    h_fluct = h_stack - mean(h_stack, 3);
    
    % 2. Calculate the variance of each pixel's fluctuations over time
    pixel_variances = var(h_fluct, 0, 3);
    
    % 3. Average the variance across the whole ROI to get total spatial-temporal power
    var_h = mean(pixel_variances(:));
    
    
    log_term = log( (sigma_mean + kappa_mean * qmax^2) / ...
        (sigma_mean + kappa_mean * qcut^2) );
    
    geometric_factor = (1 / (4 * pi)) * log_term;
    
    total_energy = sigma_mean * var_h;
    
    estimated_activity = total_energy - (kB_T)* geometric_factor;
    
    % Clean up negatives (physically impossible)
    if estimated_activity < 0; estimated_activity = 0.0032*kB_T; end
    
    map_sigma(roi_idx) = sigma_mean;
    map_activity(roi_idx) = estimated_activity/kB_T;
    
end




%%

disp('Starting ROI Analysis with 2-Parameter Fitting (New Model)...');



% Storage for the fitted parameters
fitted_eta = zeros(length(FBR), 1);
fitted_tau = zeros(length(FBR), 1);

fbr_R2_fit_d1 = zeros(length(FBR),1);
fbr_R2_fit_d2 = zeros(length(FBR),1);
fbr_R2_fit_d3 = zeros(length(FBR),1);

% Define Physics Constants
kB = 1.380649e-23;
T_kelvin = 310; % Body temp approx
kB_T = kB * T_kelvin;

% ---------------- FITTING USING IRM LOGIC ----------------

fitted_eta  = zeros(length(FBR),1);
fitted_tau  = zeros(length(FBR),1);

kB = 1.380649e-23;
T_kelvin = 310;

for roi_idx = 1: size(FBR,1)
    
    % --- ROI ---
    x = FBR(roi_idx,1);
    y = FBR(roi_idx,2);
    side_length = FBR(roi_idx,3);
    
    roi_stack = image_stack(y:y+side_length, x:x+side_length, :);
    h_stack = (double(roi_stack) - I_min) / k_factor;
    
    [Ny,Nx,Nt] = size(h_stack);
    nfft_global = 2^nextpow2(Nt);
    
    [~,f_dummy] = pwelch(squeeze(h_stack(1,1,:)), ...
        hann(floor(Nt/2)),[],nfft_global,fs);
    
    all_PSD = zeros(side_length,side_length,length(f_dummy));
    
    for k = 1:side_length
        for j = 1:side_length
            h_series = squeeze(h_stack(k,j,:));
            h_series = h_series - mean(h_series);
            
            [psd_tmp,freq_Hz_global] = pwelch(h_series,...
                hann(floor(Nt/2)),[],nfft_global,fs);
            
            all_PSD(k,j,:) = psd_tmp';
        end
    end
    
    averaged_PSD = squeeze(mean(all_PSD,[1 2],'omitnan'));
    idx_fitting = min(find(freq_Hz_global >= 0.5));
    averaged_PSD = averaged_PSD(1:idx_fitting);
    freq_Hz_global = freq_Hz_global(1:idx_fitting);   
    averaged_PSD = averaged_PSD - averaged_bg_PSD(1:idx_fitting);
    %-------- FIXED PHYSICAL PARAMETERS --------
    averaged_PSD_fitting = averaged_PSD;
    freq_Hz_global_fitting = freq_Hz_global;
    kappa_mean = 12*(kB*T_kelvin);
    sigma_mean = map_sigma(roi_idx);
    mean_confinement = 0;
    if ~isnan(sigma_mean)
        flicker_activity = map_activity(roi_idx);
        
        q = linspace(qmin,qmax,800);
        qcut = qmax;
        
        % -------- FIT PARAMETERS --------
        
        initial_guess = [1000, 30];   % [eta_eff, tau_flicker]
        lb = [1, 50e-3];
        ub = [1e6, 102.4];
        
        fit_func = @(x) IRM_cost(exp(x), ...
            freq_Hz_global_fitting, averaged_PSD_fitting, ...
            kB, T_kelvin, q, ...
            kappa_mean, sigma_mean, mean_confinement, ...
            flicker_activity, qcut, lb, ub);
        
        options = optimset('Display','off','TolX',1e-4);
        best_log = fminsearch(fit_func, log(initial_guess), options);
        best_p = exp(best_log);
        
        eta_eff     = best_p(1);
        tau_flicker = best_p(2);
        
        fitted_eta(roi_idx) = eta_eff;
        fitted_tau(roi_idx) = tau_flicker;
        
                
        % -------- GENERATE THEORETICAL PSD --------
        
        PSD_theoretical_fitting = zeros(size(freq_Hz_global_fitting));
        
        for i = 1:length(freq_Hz_global_fitting)
            
            omega = 2*pi*freq_Hz_global_fitting(i);
            
            denom_freq = ((4*eta_eff)*omega).^2;
            
            
            denom_mech_1 = ( ...
                kappa_mean .* q.^3 + ...
                sigma_mean .* q).^2;
            
            % Passive
            integrand_2 = 1 ./ (omega.^2 + denom_mech_1);
            
            % Active spatial spectrum
            flicker_Da_q = flicker_activity .* exp(-(q/qcut));
            
            % Frequency-dependent active
            S_F_q_omega = flicker_Da_q ./ (1 + (omega*tau_flicker).^2);
            
            integrand_3 = 1 ./ (denom_freq+denom_mech_1) .* S_F_q_omega;
            
            PSD_theoretical_fitting(i) = ((4*eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_2) + ...
                ((4*eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_3);
        end
        
        PSD_theoretical_fitting = PSD_theoretical_fitting * 1e18;
        
        % -------- PLOT --------
        figure(1); clf;
        loglog(freq_Hz_global, averaged_PSD,'ko'); hold on;
        loglog(freq_Hz_global_fitting, PSD_theoretical_fitting,'r-','LineWidth',2);
        legend('Data','IRM Fit');
        title(['ROI ' num2str(roi_idx)]);
        drawnow;
        
        domain_ranges = [
            0.01, 0.1;   % Domain 1: Ultra-low frequency (Active)
            0.1,  0.5;    % Domain 2: Mid-range (Transition)
            0.5,   10.0;
            ];
        
        % Log-space R^2 Function
        calc_R2 = @(obs, pred) 1 - sum((log10(obs) - log10(pred)).^2) / ...
            sum((log10(obs) - mean(log10(obs))).^2);
        
        r2_results = zeros(4, 1);
        
        % --- 3. Loop Through Domains ---
        averaged_PSD = averaged_PSD';
        PSD_theoretical_fitting = PSD_theoretical_fitting';
        
        for d = 1:4
            if d==4
                r2_results(d,1) = roi_idx;
                break;
            end
            f_start = domain_ranges(d, 1);
            f_end   = domain_ranges(d, 2);
            
            % Mask for the current domain
            mask = (freq_Hz_global_fitting >= f_start & freq_Hz_global_fitting <= f_end);
            
            if any(mask)
                r2_results(d, 1) = mean(calc_R2(averaged_PSD(mask), PSD_theoretical_fitting(mask)));
            else
                r2_results(d, :) = NaN;
            end
        end
        
        fprintf('ROI %d → eta=%.3e  tau=%.3f R2 values are %.3f and %.3f\n', ...
            roi_idx, eta_eff, tau_flicker, r2_results(1,1),r2_results(2,1));


        fbr_R2_fit_d1(roi_idx,1) = r2_results(1,1);
        fbr_R2_fit_d2(roi_idx,1) = r2_results(2,1);
        fbr_R2_fit_d3(roi_idx,1) = r2_results(3,1);
    end
end

cd(directory);
save([sprintf('Properties_all'),'.mat'],'fbr_R2_fit_d1','fbr_R2_fit_d2','fbr_R2_fit_d3','fitted_eta','fitted_tau','map_activity','map_sigma');


%% --- HELPER FUNCTION FOR OPTIMIZATION (UPDATED LOGIC) ---
function cost = IRM_cost(p, freq, psd_exp, ...
    kB, T_kelvin, q, ...
    kappa_mean, sigma_mean, mean_confinement, ...
    flicker_activity, qcut, lb, ub)

    % Bound check
    if any(p < lb) || any(p > ub)
        cost = Inf;
        return;
    end
    
    eta_eff     = p(1);
    tau_flicker = p(2);
    psd_model = zeros(size(freq));
    
    % ---- Compute model PSD ----
    for i = 1:length(freq)
        omega = 2*pi*freq(i);
        denom_freq = ((4*eta_eff)*omega).^2;
        
        denom_mech_1 = ( ...
            kappa_mean .* q.^3 + ...
            sigma_mean .* q ).^2;
        
        % Passive
        integrand_2 = 1 ./ (omega.^2 + denom_mech_1);
        
        % Active
        flicker_Da_q = flicker_activity .* exp(-(q/qcut));
        S_F_q_omega  = (flicker_Da_q) ./ (1 + (omega*tau_flicker).^2);
        
        integrand_3 = 1 ./ (denom_freq + denom_mech_1) .* S_F_q_omega;
        
        psd_model(i) = ...
            ((4*eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_2) + ...
            ((4*eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_3);
    end
    
    psd_model = psd_model * 1e18;
    
    % ---- Valid mask (ignores negative points from background subtraction) ----
    valid = psd_model > 0 & psd_exp > 0;
    if sum(valid) == 0
        cost = Inf;
        return;
    end
    
    % ---- Log residuals ----
    log_diff = log10(psd_exp(valid)) - log10(psd_model(valid));
    freq_valid = freq(valid);
    
    % ---- Frequency weighting (FIXED) ----
    weights = ones(size(freq_valid));
    weights(freq_valid >= 0.5) = 0; % Force algorithm to ignore camera noise tail
    
    % ---- Weighted cost ----
    cost = sum(weights .* (log_diff.^2));
end
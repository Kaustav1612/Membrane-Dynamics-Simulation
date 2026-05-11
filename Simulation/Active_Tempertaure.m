close all;clc;
I_min    = 1224.71868;
k_factor = 58.64994;



% ---- File Import ----
old_dir = pwd;
olddir = cd('I:\Jibitesh\IRM\110124\Control\Exported\003\set');
directory = 'K:\IRM\CytoD\C\03';
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
map_T_eff = zeros(size(FBR,1), 1);
num_rois = size(FBR,1);
roi_idxs = find(fbr_R2_fit_d1>=0.75);
for i = 1:length(roi_idxs)
    fbr_idx = roi_idxs(i);
    kappa_mean=12*kB*T_kelvin;
    sigma_mean = map_sigma(fbr_idx,1);
    eta_eff = fitted_eta(fbr_idx,1);
    x = FBR(fbr_idx,1);
    y = FBR(fbr_idx,2);
    side_length = FBR(fbr_idx,3);

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
    averaged_PSD = averaged_PSD(:);
    freq_Hz_global = freq_Hz_global(:);
    averaged_PSD = averaged_PSD - averaged_bg_PSD;

    PSD_theoretical_passive = zeros(size(freq_Hz_global));

    for f = 1:length(freq_Hz_global)

        omega = 2*pi*freq_Hz_global(f);

        denom_freq = ((4*eta_eff)*omega).^2;
        q = linspace(qmin,qmax,800);

        denom_mech_1 = ( ...
            kappa_mean .* q.^3 + ...
            sigma_mean .* q).^2;

        % Passive
        integrand_2 = 1 ./ (omega.^2 + denom_mech_1);


        PSD_theoretical_passive(f) = ((4 * eta_eff)*kB*T_kelvin/pi) * trapz(q, integrand_2);
    end

    PSD_theoretical_passive = PSD_theoretical_passive * 1e18;

    T_eff = 310*(averaged_PSD./PSD_theoretical_passive);
    idx_low_freq = min(find(freq_Hz_global>=1));
    T_eff = T_eff(1:idx_low_freq);
    avg_T_eff = mean(T_eff(:));

    map_T_eff(fbr_idx) = avg_T_eff;
end


%%
cd(old_dir);
slopes  = NaN(size(FBR,1),1);

for roi_idx = 1:size(FBR,1)
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
    results = modewise_psd_analysis(h_stack_smooth, grid_spacing, 3e7);
    slopes(roi_idx,1)=results.slope;

    q_avg = results.q;
    psd_avg = results.psd;
    q_fit = results.q_fit;
    psd_fit = results.psd_fit;
    slope = results.slope;

    if slope <= -3.5 || slope > -1

        figure;
        loglog(q_avg, psd_avg, 'bo'); hold on;
        loglog(q_fit, psd_fit, 'r-', 'LineWidth', 4);

        xlabel('q (1/m)');
        ylabel('<|h(q)|^2>');
        title(['Mode-wise PSD (Patch), slope = ' num2str(slope, '%.2f')]);
        grid on;
       
    end

end
cd(directory);
save([sprintf('active_param'),'.mat'],'slopes','map_T_eff');
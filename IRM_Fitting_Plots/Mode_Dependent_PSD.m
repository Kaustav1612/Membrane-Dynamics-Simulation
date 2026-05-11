close all;clc;
I_min    = 1382.88079;
k_factor = 46.78692;


%%
% ---- File Import ----
old_dir = pwd;
olddir = cd('/Volumes/TOSHIBA EXT/Kaustav/IRM/240126/Control/001/set');
directory = '/Volumes/TOSHIBA EXT/Kaustav/IRM/240126/Control/001/';
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


% Tension from excess area calculation for FBRs
% Tension from excess area calculation for FBRs
pixel_size = grid_spacing;
kappa_mean = 12*kB_T;

% --- 0. Setup the qcut Sweep Array ---
num_q_steps = 100;


% Create qcut sweep from low to high, avoiding qmax
qmin = 2*pi / (side_length * pixel_size); % Minimum q from system size



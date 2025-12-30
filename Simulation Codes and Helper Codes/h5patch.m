clearvars -except height_time_series
% If 'height_time_series' already exists in workspace, you can skip to section 2.
nfiles = size(height_time_series,3);

% Parameters for loading
xstart = 1; xend = 60;
ystart = 1; yend = 60;
save_dir = 'D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\protein_present\active\non_pinned\non_protrusion';
if ~exist('height_time_series', 'var')
    % Pre-allocate for speed (60x60xN)
    height_time_series = double(zeros(60, 60, nfiles));
    
    for n = 1:nfiles
        A1a = double(imread(list1(n).name));
        
        % Normalize/Crop logic from your snippet
        % (Adjusting indices to ensure we get exactly 60x60)
        A1b = (A1a(xstart:xend, ystart:yend) - inmin) ./ conv;
        
        height_time_series(:,:,n) = A1b;
    end
end



% --- 2. Extract 12x12x500 Patches ---
nn = 'Dataset_Patches'; % Output folder name
mkdir(nn);

M = 12;           % Spatial Patch Size (12x12)
T_chunk = 1000;     % Temporal Chunk Size (1000 frames)
a = height_time_series;
[H, W, T_total] = size(a);

% Create counters for loops
% 1. Spatial Row (x)
% 2. Spatial Col (y)
% 3. Time Segment (t)

patch_count = 0;

for t = 1:T_chunk:(T_total - T_chunk + 1)
    s=0;
    for y = 1:M:(W - M + 1)   
        for x = 1:M:(H - M + 1)
            s=s+1;
            patch_count = patch_count + 1;
            
            % Define ranges
            x_range = x : x + M - 1;
            y_range = y : y + M - 1;
            t_range = t : t + T_chunk - 1;
            
            % Extract the 3D block (12 x 12 x 50)
            mydata = height_time_series(x_range, y_range, t_range);
            
            % --- Save to H5 ---
            % Construct dynamic filename including spatial and time coords
            % Example: Patch_X001_Y001_T001.h5
            f = sprintf('T%04d_%d.h5',t,s);
            filename2 = fullfile(nn, f);
            
            % Create H5 file with dimensions [12 12 50]
            % Note: Check if file exists to avoid overwrite errors during testing
            if exist(filename2, 'file')
                delete(filename2);
            end
            
            h5create(filename2, '/dataset4', [M M T_chunk]);
            h5write(filename2, '/dataset4', mydata);
            
        end
    end
end

fprintf('Done. Created %d patches in folder "%s".\n', patch_count, nn);
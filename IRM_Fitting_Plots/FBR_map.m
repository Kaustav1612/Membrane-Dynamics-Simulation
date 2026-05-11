clc;


% Define range for color mapping
min_val = 1;
max_val = 2;

% Use 'jet' colormap
num_colors = 256;
c_map = jet(num_colors);

figure('Color', 'w');
% --- 3. Plot Base Image ---
imshow(A',[]);
colormap(gray);
hold on;

% --- 4. Loop through ROIs ---
num_rois = size(FBR,1);
roi_idxs = find(fbr_R2_fit_d1>=0.75);
non_fit_idx = find(fbr_R2_fit_d1>=0 & fbr_R2_fit_d2>=0.75);
for i = 1:length(roi_idxs)
    % A. Get the Value for Color (Column 2)
    fbr_idx = roi_idxs(i);
    val = map_activity(roi_idxs(i),1);
    % Map value to colormap index [1, 256]
    if max_val == min_val
        idx = 1;
    else
        % Normalize (0 to 1) -> Scale to (1 to 256)
        norm_val = (val - min_val) / (max_val - min_val);
        idx = round(norm_val * (num_colors - 1)) + 1;
    end
    
    % Clip index to be safe
    idx = max(1, min(num_colors, idx));
    
    % Get RGB color for this rectangle
    rect_color = c_map(idx, :);
    
    % B. Get Position [x, y, w, h]
    % Adjust indices based on your FBR structure:
    % Case 1: If FBR is [x, value, w, h]
    x = FBR(fbr_idx, 1);
    % y is missing? If y is separate, add it here.
    % Assuming FBR might be [x, y, w, h] but you said value is col 2?
    % If FBR is [x, y, w, h] and you want color from y, use this:
    y = FBR(fbr_idx, 2);
    w = FBR(fbr_idx, 3);
    h = FBR(fbr_idx, 4);
    
    % Standard MATLAB Rectangle Position: [x, y, w, h]
    pos = [x, y, w, h];
    
    % C. Draw Rectangle
    rectangle('Position', pos, ...
        'FaceColor', rect_color);    
end
%%
hold on;
for i = 1:length(non_fit_idx)
    % A. Get the Value for Color (Column 2)
    fbr_idx = fbrTen(non_fit_idx(i),1);
    val = map_activity(non_fit_idx(i),1);
    % Map value to colormap index [1, 256]
    if max_val == min_val
        idx = 1;
    else
        % Normalize (0 to 1) -> Scale to (1 to 256)
        norm_val = (val - min_val) / (max_val - min_val);
        idx = round(norm_val * (num_colors - 1)) + 1;
    end
    
    % Clip index to be safe
    idx = max(1, min(num_colors, idx));
    
    % Get RGB color for this rectangle
    rect_color = c_map(idx, :);
    
    % B. Get Position [x, y, w, h]
    % Adjust indices based on your FBR structure:
    % Case 1: If FBR is [x, value, w, h]
    x = FBR(fbr_idx, 1);
    % y is missing? If y is separate, add it here.
    % Assuming FBR might be [x, y, w, h] but you said value is col 2?
    % If FBR is [x, y, w, h] and you want color from y, use this:
    y = FBR(fbr_idx, 2);
    w = FBR(fbr_idx, 3);
    h = FBR(fbr_idx, 4);
    
    % Standard MATLAB Rectangle Position: [x, y, w, h]
    pos = [x, y, w, h];
    
    % C. Draw Rectangle
    rectangle('Position', pos, ...
        'FaceColor', rect_color);    
end


% --- 5. Add Colorbar ---
% This trick forces the colorbar to show the ROI scale, not the image scale
c = colorbar;
colormap(c, jet);
ylabel(c, 'x 10^3 Pa.s');

hold off;
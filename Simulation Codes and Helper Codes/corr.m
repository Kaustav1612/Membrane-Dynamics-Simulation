% === Input ===
% h: height map stack of size [Nx, Ny, Nt]

clearvars -except height_time_series
    h = height_time_series;
% refFrame: reference frame index
refFrame = 20;
[Nx, Ny, Nt] = size(h);

% --- Remove mean across time (optional but helpful) ---
h_centered = h - mean(h, 3);

% --- Reference frame ---
h_ref = h_centered(:, :, refFrame);
h_ref_temp = h_centered(:, :, 1);
%% ===============================================================
% (1) TEMPORAL AUTOCORRELATION MAP for all frames w.r.t. refFrame
% ===============================================================
temporal_corr_map = zeros(Nx, Ny, Nt);

for t = 1:Nt
    num = h_centered(:, :, t) .* h_ref;
    den = sqrt(h_centered(:, :, t).^2 .* h_ref.^2);
    temporal_corr_map(:, :, t) = num ./ (den + eps); % elementwise correlation
end

% Clip values to [-1, 1]
temporal_corr_map = max(min(temporal_corr_map, 1), -1);

% --- Example: visualize correlation map for selected frames ---
figure;
frames_to_show = [refFrame, refFrame+5, refFrame+10, Nt];
for i = 1:length(frames_to_show)
    subplot(1, length(frames_to_show), i);
    imagesc(temporal_corr_map(:, :, frames_to_show(i)), [-1 1]);
    axis image off;
    title(['Frame ' num2str(frames_to_show(i))]);
    colorbar;
end
sgtitle('Temporal autocorrelation map vs reference frame 21');


%% ===============================================================
% Average temporal correlation over all pixels for summary plot
% ===============================================================
mean_temporal_corr = squeeze(mean(mean(temporal_corr_map, 1), 2));

mean_corr_map = squeeze(mean(temporal_corr_map, 3));
figure;
imagesc(mean_corr_map);
colormap jet;
colorbar;
title('Pixel wise correlation averaged across time');

figure;
plot(1:Nt, mean_temporal_corr, '-', 'LineWidth', 1.5);
xlabel('Frame number');
ylabel('Mean temporal correlation');
title(['Average temporal correlation vs frame (ref = ' num2str(refFrame) ')']);
grid on;
%%
% Preallocate
Nt = 10000;
mean_corr = zeros(1, Nt);

% Compute temporal autocorrelation for each lag
for tau = 0:Nt
    num = 0;
    den = 0;
    count = 0;

    % Loop over all reference frames that allow lag tau
    for t = 1:(Nt - tau)
        h_t = h_centered(:, :, t);
        h_tau = h_centered(:, :, t + tau);

        % Compute dot product (sum over all pixels)
        num = num + sum(h_t(:) .* h_tau(:));
        den = den + sqrt(sum(h_t(:).^2) * sum(h_tau(:).^2));
        count = count + 1;
    end

    % Mean over all valid reference frames
    mean_corr(tau + 1) = num / (den + eps);
end

% Normalize correlation to start from 1
mean_corr = mean_corr / mean_corr(1);
%%
% Plot the result
figure;
plot(0:Nt, mean_corr, 'LineWidth', 2);
xlabel('Time lag (\tau)');
ylabel('Mean temporal autocorrelation');
title('Mean temporal autocorrelation function of height field');
grid on;

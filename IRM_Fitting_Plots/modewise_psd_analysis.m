function results = modewise_psd_patch(h, dx, max_q)

% h : (Nx x Ny x Nt) height field
% dx: spatial resolution
% max_q: maximum q to include in fit

[Nx, Ny, Nt] = size(h);

Lx = Nx * dx;
Ly = Ny * dx;

%% -------------------------------
% 1. Fourier transform
%% -------------------------------
hq = fftshift(fft2(h));  % shift zero freq to center

%% -------------------------------
% 2. Power spectrum (time averaged)
%% -------------------------------
psd2D = mean(abs(hq).^2, 3) / Nt;

%% -------------------------------
% 3. Construct qx, qy
%% -------------------------------
kx = (-Nx/2:Nx/2-1);
ky = (-Ny/2:Ny/2-1);

[qx, qy] = meshgrid(2*pi*kx/Lx, 2*pi*ky/Ly);
q_mag = sqrt(qx.^2 + qy.^2);

%% -------------------------------
% 4. Radial averaging (binning)
%% -------------------------------
q_flat = q_mag(:);
psd_flat = psd2D(:);

% Remove zero mode
valid = q_flat > 0;
q_flat = q_flat(valid);
psd_flat = psd_flat(valid);

% Define bins
num_bins = round(sqrt(Nx*Ny));  % reasonable choice
q_bins = linspace(min(q_flat), max(q_flat), num_bins);

q_avg = zeros(num_bins-1,1);
psd_avg = zeros(num_bins-1,1);

for i = 1:num_bins-1
    mask = (q_flat >= q_bins(i) & q_flat < q_bins(i+1));
    
    if any(mask)
        q_avg(i) = mean(q_flat(mask));
        psd_avg(i) = mean(psd_flat(mask));
    else
        q_avg(i) = NaN;
        psd_avg(i) = NaN;
    end
end

% Remove NaNs
valid = ~isnan(q_avg) & ~isnan(psd_avg);
q_avg = q_avg(valid);
psd_avg = psd_avg(valid);

%% -------------------------------
% 5. Limit q range
%% -------------------------------
mask = q_avg <= max_q;
q_fit = q_avg(mask);
psd_fit_data = psd_avg(mask);

%% -------------------------------
% 6. Log-log fit
%% -------------------------------
log_q = log10(q_fit);
log_psd = log10(psd_fit_data);

p = polyfit(log_q, log_psd, 1);

slope = p(1);
psd_fit = 10.^(polyval(p, log_q));

%% -------------------------------
% 7. Store results
%% -------------------------------
results.q = q_avg;
results.psd = psd_avg;
results.q_fit = q_fit;
results.psd_fit_data = psd_fit_data;
results.psd_fit = psd_fit;
results.slope = slope;
results.fit_coeff = p;

%% -------------------------------
% 8. Plot
%% -------------------------------
figure;
loglog(q_avg, psd_avg, 'bo'); hold on;
loglog(q_fit, psd_fit, 'r-', 'LineWidth', 4);

xlabel('q (1/m)');
ylabel('<|h(q)|^2>');
title(['Mode-wise PSD (Patch), slope = ' num2str(slope, '%.2f')]);
grid on;
close;
end
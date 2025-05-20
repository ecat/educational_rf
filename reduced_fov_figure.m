addPaths()
rf_pulse_path = './pulses/';

rf_pulse_tag_helper = @(x) 'msinc_tbw_2';

rf_pulse_tag = rf_pulse_tag_helper();
[Nt, dT, ~, ~, pulse_shape_complex, max_B1_G] = ...
    load_designed_pulse(rf_pulse_path, rf_pulse_tag, 0);

gradient_rewind_scale = 2;
Nt_final = Nt * 40;
time_scale = .5; % squish pulse to .5 ms
dT = dT * time_scale;
max_B1_G = max_B1_G / time_scale / 2;

gy_1 = ones(Nt, 1); 
gz = zeros(Nt, 1);
gy = gy_1;
N_rewind_subpulse = Nt / gradient_rewind_scale;

N_subpulses = 18;
weights = pulse_shape_complex(1:round(Nt / N_subpulses):end);

triangle_shape = triang(N_rewind_subpulse);
triangle_shape = triangle_shape / max(abs(triangle_shape(:)));

B1_base = pulse_shape_complex / time_scale / sum(weights) * 4;
B1_t = B1_base * weights(1);
for ii = 2:N_subpulses
    B1_t = cat(1, B1_t, zeros(N_rewind_subpulse, 1), B1_base * weights(ii));
    gy = cat(1, gy, -gradient_rewind_scale * gy_1(1:N_rewind_subpulse), gy_1);
    gz = cat(1, gz, .9 * triangle_shape, 0 * gy_1);
end
gy(end + 1) = 0;
B1_t(end + 1) = 0;

%%

my_fig = figure('Color', 'white', 'Position', [100 100 800 400]); 
plot([0; B1_t]); 
axis off;
xlim([0 inf])
ylim([-max_B1_G * .8, max_B1_G * 1.2]) % so that different pulses are scaled p
hold on;
plot([0; gy] * 0.02, 'k-');
plot([0; gz] * 0.01, 'r-');
legendflex({'$B_1(t)$', '$G_y$', '$G_z$'}, 'nrow', 1, 'FontSize', 28, 'anchor', {'s', 's'})
%legend('$B_1(t)$', '$G_y$', '$G_z$', 'FontSize', 16)

export_fig(my_fig, sprintf('out2/reduced_fov.png'))

%%
ETL = 8;
ky_t = [];
kz_t = [];
N = 128;
for ee = 1:8
    ky_t = cat(1, ky_t, (-1) ^ ee .* vec(linspace(-1, 1, N)));
    kz_t = cat(1, kz_t, ee - 4.5 * ones(128, 1));
end


figure('Color', 'white', 'Position', [100 100 400 400])
plot(-ky_t, kz_t, 'k', 'LineWidth', 3);
axis off;
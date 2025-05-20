pulse_option = 3;

if(pulse_option == 1)
    N_subpulses = 2;
    gradient_rewind_scale = 1;
elseif(pulse_option == 2)
    N_subpulses = 5;
    gradient_rewind_scale = 1; % longest pulse
elseif(pulse_option == 3)
    N_subpulses = 5;
    gradient_rewind_scale = 2;
end

rf_pulse_path = './pulses/';

rf_pulse_tag_helper = @(x) 'msinc_tbw_2';

rf_pulse_tag = rf_pulse_tag_helper();
[Nt, dT, ~, ~, pulse_shape_complex, max_B1_G] = ...
    load_designed_pulse(rf_pulse_path, rf_pulse_tag, 0);

Nt_final = Nt * 10;
time_scale = .5; % squish pulse to .5 ms
dT = dT * time_scale;
max_B1_G = max_B1_G / time_scale / 2;

gz_1 = ones(Nt, 1); 
gz = gz_1;
N_rewind_subpulse = Nt / gradient_rewind_scale;
pascal_weights = get_pascal_vector(N_subpulses);
B1_base = pulse_shape_complex / time_scale / sum(pascal_weights);
B1_t = B1_base * pascal_weights(1);
for ii = 2:N_subpulses
    B1_t = cat(1, B1_t, zeros(N_rewind_subpulse, 1), B1_base * pascal_weights(ii));
    gz = cat(1, gz, -gradient_rewind_scale * gz_1(1:N_rewind_subpulse), gz_1);
end

B1_t = cat(1, B1_t, zeros(Nt_final - numel(B1_t), 1));
gz = cat(1, gz, zeros(Nt_final - numel(gz), 1));

dfs_to_simulate = linspace(-8, 8, 200);

N_bulk_off_resonances = 100;
bulk_off_resonances = linspace(-1.2, 1.2, N_bulk_off_resonances);

M_transverse_zf = zeros(numel(dfs_to_simulate), N_bulk_off_resonances);
parfor ff = 1:N_bulk_off_resonances
    df = bulk_off_resonances(ff);
    [M_result, G] = simulate_rfpulse_every_timepoint(dfs_to_simulate, B1_t, ...
        numel(B1_t), dT, 0, 0, gz, df);

    M_xy = M_result(1, :, end) + 1i * M_result(2, :, end);
    M_transverse_zf(:, ff) = M_xy(:);
end

%%
my_fig = figure('Color', 'white', 'Position', [100 100 800 200]); 
sp1 = subplot(121);
original_position = sp1.Position;
sp1.Position(1) = sp1.Position(1) + .13;
sp1.Position(3) = sp1.Position(3) - .1;
plot([0; B1_t]); 
axis off;
xlim([0 Nt_final])
ylim([-max_B1_G * 1.2, max_B1_G * 1.2]) % so that different pulses are scaled p
hold on;
plot([0; gz] * 0.02, 'k-');

sp2 = subplot(122);
imagesc(abs(M_transverse_zf), 'YData', dfs_to_simulate, 'XData', bulk_off_resonances);

ax_fontsize = 12;
ax = gca;
ax.XAxis.FontSize = ax_fontsize;
ax.YAxis.FontSize = ax_fontsize;
xlabel('Off-Resonance [kHz]', 'FontSize', 16);
ylabel('z [cm]', 'FontSize', 16);
title('Mxy Profile', 'Fontsize', 16);

sp1.Position(2) = sp2.Position(2) + .03;
sp1.Position(4) = sp2.Position(4);

export_fig(my_fig, sprintf('out2/spsp_%d.png', pulse_option))
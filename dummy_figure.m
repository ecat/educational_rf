rf_pulse_path = './pulses/';
dfs_to_simulate = linspace(-3, 3, 800);


rf_pulse_tag = 'msinc_tbw_8_fa_45';
[Nt, dT, ~, ~, pulse_shape_complex, max_B1_G] = ...
    load_designed_pulse(rf_pulse_path, rf_pulse_tag, 0);

B1_t = pulse_shape_complex;

pulse_duration = Nt * dT;
[M_result_1, G] = simulate_rfpulse_every_timepoint(dfs_to_simulate, B1_t, ...
    numel(B1_t), dT, pulse_duration * .51, 1, [], 0);


[M_result_2, G] = simulate_rfpulse_every_timepoint(dfs_to_simulate, B1_t, ...
    numel(B1_t), dT * 2, 2 * pulse_duration * .51, 1, [], 0);

%%
plot_option = 2;

fig1a = figure('Position', [100 100 1200 400], 'Color', 'white'); 
has = tight_subplot(1, 2, [.1 .025], [.1 .1], [.1 .1]);

colora = 'r'; % [0 .447 .741]
colorb = 'b'; % [.85 .325 0.098]

axes(has(1));
if(plot_option == 1)
    plot(linspace(0, pulse_duration, Nt), B1_t, 'Color', colora); hold on;
    plot(linspace(0, pulse_duration / 2, Nt) + pulse_duration / 4, B1_t, 'Color', colorb)
    text(pulse_duration / 2, -0.02, 'Excitation k-Space', 'FontSize', 24, 'HorizontalAlignment', 'center');
else
    plot(linspace(0, pulse_duration, Nt), B1_t, 'Color', 'k'); hold on;
    text(pulse_duration / 2, -0.02, 'Pulse', 'FontSize', 24, 'HorizontalAlignment', 'center');
    text(pulse_duration * .5, 1.15 * max_B1_G, '$B_1(t)$', 'FontSize', 24, 'HorizontalAlignment', 'center');
end
xlabel('Time [ms]', 'FontSize', 16)
ylim([-0.02, 0.065])
axis off;


axes(has(2));
if(plot_option == 1)
    plot(dfs_to_simulate, M_result_1(1, :, end), 'Color', colora, 'LineStyle', ':');  hold on;
    plot(dfs_to_simulate, M_result_1(2, :, end), 'Color', colora, 'LineStyle', '-');
    plot(dfs_to_simulate, M_result_2(1, :, end), 'Color', colorb, 'LineStyle', ':');  hold on;
    plot(dfs_to_simulate, M_result_2(2, :, end), 'Color', colorb, 'LineStyle', '-');
    text(0, -0.25, 'Slice Location', 'FontSize', 24, 'HorizontalAlignment', 'center');
else
    plot(dfs_to_simulate, M_result_1(1, :, end), 'Color', 'k', 'LineStyle', ':');  hold on;
    plot(dfs_to_simulate, M_result_1(2, :, end), 'Color', 'k', 'LineStyle', '-');
    text(0, -0.25, 'Magnetization Profile', 'FontSize', 24, 'HorizontalAlignment', 'center');
    text(0, .85, '$m(z)$', 'FontSize', 24, 'HorizontalAlignment', 'center');
end
ylim([-0.25, 1.1])
xlim([-3 3])
axis off;

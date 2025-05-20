rf_pulse_path = './pulses/';

example_index = 5;

if(example_index == 0)
    tags_to_compare = {'msinc_tbw_2', 'msinc_tbw_4', 'msinc_tbw_8'};
    titlesa = {'Windowed Sinc'};
    titlesb = {'Transverse mx / my Profiles'};
    do_mxy_together = 2;
    legend_labels = {'TBW 2', 'TBW 4', 'TBW 8'};
    dfs_to_simulate = linspace(-4, 4, 1000);
    pulse_ylims = [-inf 0.15];
elseif(example_index == 1)
    tags_to_compare = {'hard_max_b1_0.03', 'hard_max_b1_0.06', 'hard_max_b1_0.12'};
    titlesa = {'Hard Pulse 1 ms'};
    titlesb = {'Longitudinal Mz Profiles'};    
    do_mxy_together = 0;
    legend_labels = {'FA 45$^\circ$', 'FA 90$^\circ$', 'FA 180$^\circ$'};
    dfs_to_simulate = linspace(-2, 2, 1000);
    pulse_ylims = [-0.002, 0.16];
elseif(example_index == 2)
    tags_to_compare = {'wurst_b1_0.10', 'wurst_b1_0.20', 'wurst_b1_0.40'};
    legend_labels = {'B1 0.1 G', 'B1 0.2 G', 'B1 0.4 G'};
    titlesa = {'WURST Pulse'};
    titlesb = {'Longitudinal Mz Profiles'};
    do_mxy_together = 0;
    dfs_to_simulate = linspace(-14, 14, 1000);
    pulse_ylims = [-0.15, 0.15];
elseif(example_index == 3)
    tags_to_compare = {'msinc_tbw_8_fa_45', 'msinc_tbw_8_fa_90', 'msinc_tbw_8_fa_180'};
    titlesa = {'Windowed Sinc'};
    titlesb = {'Transverse mx / my Profiles'};
    do_mxy_together = 2;
    legend_labels = {'FA 45$^\circ$', 'FA 90$^\circ$', 'FA 180$^\circ$'};
    dfs_to_simulate = linspace(-4, 4, 1000);
    pulse_ylims = [ -0.0500    0.2500];
elseif(example_index == 39)
    % for animating 3
    tags_to_compare = {'msinc_tbw_8_fa_45'};
    titlesa = {'Windowed Sinc'};
    titlesb = {'Transverse mx / my Profiles'};
    do_mxy_together = 2;
    legend_labels = {'FA 45$^\circ$'};
    dfs_to_simulate = linspace(-4, 4, 1000);
    pulse_ylims = [ -0.0500    0.2500];
elseif(example_index == 4)
    tags_to_compare = {'multiband_1', 'multiband_2', 'multiband_4'};
    titlesa = {'Modulated Windowed Sinc'};
    titlesb = {'Tranverse mx / my Profiles'};
    do_mxy_together = 2;
    legend_labels = {'Singleband', 'Multiband-2', 'Multiband-4'};
    dfs_to_simulate = linspace(-8, 8, 1000);
    pulse_ylims = [-.16 .16];
elseif(example_index == 5)
    tags_to_compare = {'msinc_tbw_8_fa_180', 'sbse_slfrank'};
    titlesa = {'Excitation k-Space vs SLR Design'};
    titlesb = {'Mz Profiles'};
    do_mxy_together = 0;
    legend_labels = {'Sinc', 'SLR'};
    dfs_to_simulate = linspace(-2, 2, 1000);
    pulse_ylims = [-.15 .4];
end

for tt = 1:numel(tags_to_compare)
    rf_pulse_tag = tags_to_compare{tt};
    
    [Nt, dT, ~, ~, pulse_shape_complex, max_B1_G] = ...
        load_designed_pulse(rf_pulse_path, rf_pulse_tag, 0);


    if(strcmp(tags_to_compare{tt}, 'sbse_slfrank'))
        target_Nt = 1000;
        dT = dT * 2;
        max_B1_G = max_B1_G / 2 * target_Nt / Nt * 1.5;
        pulse_shape_complex = interp1(linspace(0, 1, Nt), pulse_shape_complex, linspace(0, 1, target_Nt)).';
        pulse_shape_complex = real(pulse_shape_complex) / max(abs(pulse_shape_complex)) * max_B1_G;
        Nt = target_Nt;
    end
    B1_t = pulse_shape_complex;

    pulse_duration = Nt * dT;
    [M_result, G] = simulate_rfpulse_every_timepoint(dfs_to_simulate, B1_t, ...
        numel(B1_t), dT, pulse_duration * .51, 1, [], 0);

    dTs{tt} = dT;
    M_results{tt} = M_result;
    B1s{tt} = B1_t;
    Gs{tt} = G;
    max_B1_Gs{tt} = max_B1_G;
    N_timepoints{tt} = Nt;
end

%%
fig = figure('Position', [100 100 1200 400], 'Color', 'white');

% pad so that the plots look nicer
N_pad_left = round(Nt / 20);
N_pad_right = round(Nt / 20);

ax_fontsize = 14;
if(example_index == 5)
    linestyles = {'k-', 'r-'};
else
    linestyles = {'r-', 'b-', 'k-'};
end

pad_helper1 = @(x) padarray(x, N_pad_left, 0, 'pre');
pad_helper2 = @(x) padarray(x, N_pad_right, 0, 'post');
pad_helper = @(x) pad_helper2(pad_helper1(x));

subplot(121);

if(example_index == 2)
    ts = linspace(-dT * N_pad_left, dT * (Nt + N_pad_right), Nt + (N_pad_left + N_pad_right));
    plot(ts, abs(pad_helper(B1s{1})), 'k-', 'LineWidth', 2.5);
    yyaxis right;
    plot(ts, angle(pad_helper(B1s{1})), ':', 'LineWidth', 1.5);
    ylabel('Phase [rad]', 'FontSize', 18);
    ylim([-pi, pi])

    ax = gca;
    ax.XAxis(1).FontSize = ax_fontsize;
    ax.YAxis(1).FontSize = ax_fontsize;

    yyaxis left
    ax.XAxis(1).FontSize = ax_fontsize;
    ax.YAxis(1).FontSize = ax_fontsize;
else
    for tt = 1:numel(tags_to_compare)
        ts = linspace(-dTs{tt} * N_pad_left, dTs{tt} * (Nt + N_pad_right), Nt + (N_pad_left + N_pad_right));
        plot(ts, pad_helper(B1s{tt}), ...
            linestyles{tt}, 'LineWidth', 2.5); hold on;
    end
    ax = gca;
    ax.XAxis.FontSize = ax_fontsize;
    ax.YAxis.FontSize = ax_fontsize;
end

if(example_index == 5)
    legend(legend_labels, 'Location', 'south east', 'FOntSize', 14)
end
ylim(pulse_ylims)
xlim([min(ts) max(ts)])
ylabel('B1 [G]', 'FontSize', 18)
xlabel('Time [ms]', 'FontSize', 18)
title(titlesa, 'FontSize', 20)

if(do_mxy_together == 2)
    if(example_index == 4)
        legend(legend_labels, 'FontSize', 14, 'Location', 'south west', 'Interpreter', 'latex')
    else
        legend(legend_labels, 'FontSize', 14, 'Location', 'north east', 'Interpreter', 'latex')
    end
end
subplot(122);

for tt = 1:numel(tags_to_compare)
    if(do_mxy_together == 1)
        plot(dfs_to_simulate, sos(M_results{tt}(1:2, :, end), 1).', linestyles{tt}, 'LineWidth', 2.5); hold on;
    elseif(do_mxy_together == 2)
        plot(dfs_to_simulate, M_results{tt}(2, :, end), linestyles{tt}, 'LineWidth', 2); hold on;
        plot(dfs_to_simulate, M_results{tt}(1, :, end), strcat(linestyles{tt}(1), ':'), 'LineWidth', 1.5); 
    else
        plot(dfs_to_simulate, M_results{tt}(3, :, end).', linestyles{tt}, 'LineWidth', 2.5); hold on;

    end
end

ax = gca;
ax.XAxis.FontSize = ax_fontsize;
ax.YAxis.FontSize = ax_fontsize;
ylabel({'Magnetization'}, 'FontSize', 18)    
if(do_mxy_together == 2)
    legend('mx', 'my', 'FontSize', 14, 'Interpreter', 'latex')
else
    legend(legend_labels, 'FontSize', 14, 'Location', 'south east', 'Interpreter', 'latex')
end
ylim([-1 1])
xlim([-inf inf])
xlabel('Frequency [kHz]', 'FontSize', 18)
title(titlesb, 'FontSize', 20)


export_fig(fig, sprintf('out2/static_pulse_profiles_example_%d.png', example_index), '-nocrop', '-m2.5');

%%
if(example_index == 2)
    figb = figure('Position', [900 100 600 400], 'Color', 'white');
    imagesc(squeeze(M_results{2}(3, :, :)), 'XData', ts,'YData', dfs_to_simulate); 
    cbar = colorbar;
    cbar.TickLabelInterpreter = 'latex';
    cbar.FontSize = 16;
    %ylabel(cbar, 'Mz', 'FontSize', 22, 'Interpreter', 'latex')
    caxis([-1, 1]);
    ax = gca;
    ax.XAxis.FontSize = ax_fontsize;
    ax.YAxis.FontSize = ax_fontsize;
    xlabel('Time [ms]', 'FontSize', 18);
    ylabel('Frequency [kHz]', 'FontSize', 18);
    title('Mz During Frequency Sweep', 'FontSize', 20)
    export_fig(figb, sprintf('out2/static_pulse_profiles_example_%d_time.png', example_index), '-nocrop');
end
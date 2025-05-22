addPaths()

%%
N_mb_indices = 1;
M_results = cell(N_mb_indices, 1);
B1s = cell(N_mb_indices, 1);
Gs = cell(N_mb_indices, 1);
dTs = cell(N_mb_indices, 1);
max_B1_Gs = cell(N_mb_indices, 1);
N_timepoints = cell(N_mb_indices, 1);

pulse_option = 14;
line_to_highlight = 'mxy';
dfs_to_simulate = linspace(-5, 5, 500);
keypoints = [];
if(N_mb_indices == 1)
    rf_pulse_path = './pulses/';


    if(pulse_option == 0)
        rf_pulse_tag_helper = @(x) 'sbse_slfrank';
        line_to_highlight = 'mz';
    elseif(pulse_option == 1)
    elseif(pulse_option == 2)
        rf_pulse_tag_helper = @(x) 'hard_max_b1_0.03';
    elseif(pulse_option == 3)
        rf_pulse_tag_helper = @(x) 'hard_max_b1_0.06';
    elseif(pulse_option == 4)
        rf_pulse_tag_helper = @(x) 'hard_max_b1_0.12';
    elseif(pulse_option == 5)
        rf_pulse_tag_helper = @(x) 'msinc_tbw_2';
    elseif(pulse_option == 6)
        rf_pulse_tag_helper = @(x) 'msinc_tbw_4';
    elseif(pulse_option == 61)
        rf_pulse_tag_helper = @(x) 'msinc_truncated_tbw_4';
    elseif(pulse_option == 7)
        rf_pulse_tag_helper = @(x) 'msinc_tbw_8';
    elseif(pulse_option == 8)
        rf_pulse_tag_helper = @(x) 'msinc_tbw_8_fa_45';
    elseif(pulse_option == 9)
        rf_pulse_tag_helper = @(x) 'msinc_tbw_8_fa_90';
    elseif(pulse_option == 10)
        rf_pulse_tag_helper = @(x) 'msinc_tbw_8_fa_180';
        line_to_highlight = 'mz';
    elseif(pulse_option == 11)
        rf_pulse_tag_helper = @(x) 'wurst_b1_0.10';
        dfs_to_simulate = linspace(-10, 10, 500);
        line_to_highlight = 'mz';
    elseif(pulse_option == 12)
        rf_pulse_tag_helper = @(x) 'wurst_b1_0.20';
        dfs_to_simulate = linspace(-10, 10, 500);
        line_to_highlight = 'mz';
    elseif(pulse_option == 13)
        rf_pulse_tag_helper = @(x) 'wurst_b1_0.40';
        dfs_to_simulate = linspace(-10, 10, 500);
        line_to_highlight = 'mz';
    elseif(pulse_option == 14 || pulse_option == 15) % binomial pulse examples
        rf_pulse_tag_helper = @(x) 'msinc_tbw_8_fa_45';
    elseif(pulse_option == 16)
        rf_pulse_tag_helper = @(x) 'multiband_1';
        dfs_to_simulate = linspace(-15, 15, 500);
    elseif(pulse_option == 17)
        rf_pulse_tag_helper = @(x) 'multiband_2';
        dfs_to_simulate = linspace(-15, 15, 500);
    elseif(pulse_option == 18)
        rf_pulse_tag_helper = @(x) 'multiband_4';
        dfs_to_simulate = linspace(-15, 15, 500);
    else
        error()
    end
end


for mb_index = 1:N_mb_indices   

    rf_pulse_tag = feval(rf_pulse_tag_helper, mb_index);
    do_plot_pulse = 0;
    [Nt, dT, ~, ~, pulse_shape_complex, max_B1_G] = ...
        load_designed_pulse(rf_pulse_path, rf_pulse_tag, do_plot_pulse);

    
    iso2end_time = dT * Nt / 2 * 1.01;
    N_rewinder = round(Nt * .45);
    G = [];
    constant_off_resonance = 0;
    if(pulse_option == 0)
        dT_new = 4 / Nt;
        max_B1_G = max_B1_G * dT / dT_new;
        pulse_shape_complex = pulse_shape_complex * dT / dT_new;
        dT = dT_new;
    elseif(pulse_option == 6)
        keypoints = ones(Nt, 1);
        keypoints = cat(1, keypoints, 2 * ones(N_rewinder, 1));
        %iso2end_time = dT * Nt * 0.55; % manually tuned
        iso2end_time = dT * Nt * .5;
        %pulse_shape_complex = pulse_shape_complex * .5;
        %max_B1_G = max_B1_G * .5;
    elseif(pulse_option == 61)
        N_rewinder = round(.45 * Nt);
        %iso2end_time = dT * Nt * .28; manually tuned
        iso2end_time = dT * Nt * .25;
        %pulse_shape_complex = pulse_shape_complex * .5;
        %max_B1_G = max_B1_G * .5;
        keypoints = ones(Nt, 1);
        keypoints = cat(1, keypoints, 2 * ones(N_rewinder, 1));
    elseif(pulse_option == 14 || pulse_option == 15)
        pulse_shape_complex = cat(1, pulse_shape_complex, zeros(Nt, 1), pulse_shape_complex, zeros(Nt/2, 1));
        G = cat(1, ones(Nt, 1), -ones(Nt, 1), ones(Nt, 1), -ones(Nt/2, 1));
        iso2end_time = 0;
        N_rewinder = 0;

        if(pulse_option == 15)
            peak_to_peak_separation_ms = 8;
            constant_off_resonance = 1/(2 * 8);
        end

        keypoints = ones(2 * Nt, 1);  % frame indices to start a new file
        keypoints = cat(1,keypoints, 2 * ones(numel(pulse_shape_complex) - numel(keypoints), 1));
    end

    B1_t = pulse_shape_complex;
    [M_result, G] = simulate_rfpulse_every_timepoint(dfs_to_simulate, B1_t, ...
        numel(B1_t), dT, iso2end_time, N_rewinder, G, constant_off_resonance);
    dTs{mb_index} = dT;
    M_results{mb_index} = M_result;
    B1s{mb_index} = B1_t;
    Gs{mb_index} = G;
    max_B1_Gs{mb_index} = max_B1_G;
    N_timepoints{mb_index} = Nt;
end

dT = dTs{1};
max_B1_G = max_B1_Gs{1};
N_timepoints = N_timepoints{1};
%%
for mb_index = 1:N_mb_indices
    M_result = M_results{mb_index};
    M_xy_t = squeeze(M_result(1, :, :) + 1i * M_result(2, :, :));

    figure; 
    subplot(311);
    imagesc(abs(M_xy_t)); colormap gray;
    xlabel('t')
    ylabel('f [kHz]');
    title('|Mxy|')

    subplot(312);
    imagesc(angle(M_xy_t)); colormap gray;
    xlabel('t')
    ylabel('f [kHz]');
    title('angle Mxy')
    
    subplot(313);
    plot(dfs_to_simulate, abs(M_xy_t(:, end)));
    yyaxis right;
    plot(dfs_to_simulate, angle(M_xy_t(:, end)));
    xlabel('f [kHz]')
    
    suptitle(sprintf('Pulse %d', mb_index - 1))

end

%%
mb_index_to_plot = 1;
do_anim = 1;
do_slice_select = ~(any(pulse_option == [2 3 4 12]));
do_overlay_signal = do_slice_select && any(pulse_option == [5 6 7 61]);
do_show_all_m = any(pulse_option == [1 5 6 7 14 15 16 17 61]);
do_show_signal_title = 1 && any(pulse_option == [1 6 61]);

G_amplitude = 1; 


normalization_factor = abs(sum(M_result(1, :, end) + 1i * M_result(2, :, end))); % so that last is 1.0
measured_signal_as_function_of_t = squeeze(sum(M_result(1, :, :) + 1i * M_result(2, :, :), 2) / normalization_factor);

if(do_slice_select)
    N_total_points = size(M_result, 3);
    G_to_plot = Gs{mb_index_to_plot};
    B1_to_plot = B1s{mb_index_to_plot};
    pulse_profile_xvals = dfs_to_simulate / (4.258 * G_amplitude);
    plot1_xlabel = 'Location $z = f / (G_z\gamma)$ [cm]';
    T_total = dT * N_total_points;
    zpad_helper = @(x) padarray(x, N_total_points - numel(x), 0, 'post');
    ts = linspace(0, T_total, N_total_points);
else
    N_total_points = N_timepoints;
    G_to_plot = zeros(N_total_points, 1);
    B1_to_plot = B1s{mb_index_to_plot}(1:N_total_points);
    pulse_profile_xvals = dfs_to_simulate;
    plot1_xlabel = 'Frequency $f$ [kHz]';
    plot1_xticks = [-5, -2.5, 0, 2.5, 5];
    T_total = dT * N_total_points;
    zpad_helper = @(x) x;
    ts = linspace(0, T_total, N_total_points);
end

if(pulse_option == 2 || pulse_option == 3 || pulse_option == 4)
    timepoints_to_animate = 1:1:N_total_points;
elseif(pulse_option == 6)
    timepoints_to_animate = 1:10:N_total_points;
elseif(pulse_option == 61)
    timepoints_to_animate = 1:10:N_total_points;
elseif(pulse_option == 12)
    timepoints_to_animate = 1:10:N_total_points;
elseif(pulse_option == 14 || pulse_option == 15)
    timepoints_to_animate = 1:30:N_total_points;
else
    timepoints_to_animate = 1:30:N_total_points;
end

highlight_color = [.635 .079 .184];

if(strcmp(line_to_highlight, 'mz'))
    mz_color = highlight_color;
    mxy_color = 'k';
else
    mz_color = 'k';
    mxy_color = highlight_color;
end

M_result = M_results{mb_index_to_plot};
fig = figure('Position', [100 100 900 500], 'Color', 'white');

sp1 = subplot(1, 2, 1);
axis([pulse_profile_xvals(1) pulse_profile_xvals(end) -1.01 1.01]);
h2b = animatedline('MaximumNumPoints', numel(dfs_to_simulate), 'LineWidth', 2.5, 'Color', mxy_color, 'LineStyle', '-');
if(do_show_all_m)
    h1 = animatedline('MaximumNumPoints', numel(dfs_to_simulate), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
    h2 = animatedline('MaximumNumPoints', numel(dfs_to_simulate), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
end
h3 = animatedline('MaximumNumPoints', numel(dfs_to_simulate), 'LineWidth', 2.5, 'Color', mz_color, 'LineStyle', ':');
ax_fontsize = 14;
label_fontsize = 18;
if(do_show_all_m)
    legend('$|$Mxy$|$', 'Mx', 'My', 'Mz', 'Location', 'south west', 'FontSize', 12);
else
    legend('$|$Mxy$|$', 'Mz', 'Location', 'south west', 'FontSize', 12);
end
ax = gca;
ax.XAxis.FontSize = ax_fontsize;
ax.YAxis.FontSize = ax_fontsize;
ylabel('Magnetization [a.u.]', 'FontSize', label_fontsize);
xlabel(plot1_xlabel, 'FontSize', label_fontsize);

if(do_show_signal_title)
    txt_handle = text(0, 1.05, signal_label_getter(0), 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

plot_zero_line = @() plot(ts, zeros(size(ts)), 'Color', [.9 .9 .9]);

sp2 = subplot(2, 2, 2);
axis([0 T_total -max_B1_G max_B1_G])
plot(ts, zpad_helper(real(B1_to_plot)), 'k-'); hold on;
plot(ts, zpad_helper(imag(B1_to_plot)), 'k:'); 
xlim([0, T_total])
if(pulse_option == 12)
    ylim([-inf inf]) 
else
    ylim([-.5 * max_B1_G 1.03 * max_B1_G]) 
end
%lgd = legend('Real', 'Imag', 'AutoUpdate', 'off', 'FontSize', 12);
%lgd.Position(1) = lgd.Position(1) - .01;
plot_zero_line();
set(gca, 'Children', flipud(get(gca, 'Children')) )
h4 = animatedline('MaximumNumPoints', 2, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
ax = gca;
ax.XAxis.FontSize = ax_fontsize;
ax.YAxis.FontSize = ax_fontsize;
ylabel('Pulse $B_1(t)$ [G]', 'FontSize', label_fontsize);

if(do_overlay_signal)
    yyaxis right;
    h_signal_overlay = animatedline('MaximumNumPoints', N_total_points, 'Color', [.494 .184 .556] + .1, 'LineWidth', 2.5);
    %plot(ts(1:N_total_points), imag(measured_signal_as_function_of_t(1:N_total_points)), 'Color', [0.494 .184 .556]);
    ylabel('Total Signal', 'FontSize', label_fontsize)
    ax = gca;
    ax.YAxis(2).FontSize = ax_fontsize;
    % keep same aspect ratio as pulse
    ylim([-.5, 1.01])
    yyaxis left;
end

sp3 = subplot(2, 2, 4);
if(pulse_option == 6 || pulse_option == 61)
    axis([0, round(T_total), -1.3 1.3]);
else
    axis([0, round(T_total), -max(abs(Gs{1})) max(abs(Gs{1}))]);
end
plot(ts, zpad_helper(G_to_plot), 'k'); hold on;
plot_zero_line();
set(gca, 'Children', flipud(get(gca, 'Children')) )
h5 = animatedline('MaximumNumPoints', 2, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
if(pulse_option == 6 || pulse_option == 61)
    ylim([-1.3, 1.3])
    xlim([0, T_total])
else
    ylim(1.1 * [-max(abs(Gs{1})) max(abs(Gs{1}))])
    xlim([0, T_total])
end
ax = gca;
ax.XAxis.FontSize = ax_fontsize;
ax.YAxis.FontSize = ax_fontsize;
ylabel('$G_z$ [G/cm]', 'FontSize', label_fontsize)
xlabel('Time [ms]', 'FontSize', label_fontsize)

%
do_save_fig = 2;
curr_keypoint_index = 0;

for t = 1:numel(timepoints_to_animate)

    tt = timepoints_to_animate(t);

    if(numel(keypoints) == 0)
        keypoint_index = 1;
    else
        keypoint_index = keypoints(tt);
    end
    if(keypoint_index > curr_keypoint_index)
        if(curr_keypoint_index > 0)
           close(v)
        end
        v = VideoWriter(sprintf('out2/profile_evolution_%d_step_%d.avi', pulse_option, curr_keypoint_index), 'Motion JPEG AVI');
        v.Quality = 75;

        if(pulse_option == 14 || pulse_option == 15)
            v.FrameRate = 15;
        end
        open(v)
        curr_keypoint_index = keypoint_index;
    end

    if(do_show_all_m)
        clearpoints(h1)
        clearpoints(h2)
        addpoints(h1, pulse_profile_xvals, squeeze(M_result(1, :, tt)))
        addpoints(h2, pulse_profile_xvals, squeeze(M_result(2, :, tt)));
    end
    addpoints(h2b, pulse_profile_xvals, abs(M_result(1, :, tt) + 1i * M_result(2, :, tt)));
    clearpoints(h3)
    addpoints(h3, pulse_profile_xvals, squeeze(M_result(3, :, tt)));

    if(do_show_signal_title)
        txt_handle.String = signal_label_getter(measured_signal_as_function_of_t(tt));
    end

    if(do_overlay_signal)
        clearpoints(h_signal_overlay);
        addpoints(h_signal_overlay, ts(1:tt), imag(measured_signal_as_function_of_t(1:tt)));
    end

    clearpoints(h4)
    addpoints(h4, ones(2, 1) * ts(tt), sp2.YLim)

    clearpoints(h5)
    addpoints(h5, ones(2, 1) * ts(tt), sp3.YLim)

    if(do_save_fig == 1)
        drawnow;
        out_folder = sprintf('./out2/%s_ss_%d_overlay_signal_%d/', rf_pulse_tag_helper(mb_index_to_plot), do_slice_select, do_overlay_signal);
        out_fig_path = sprintf('%s/frame_%d.png', out_folder, t);
    
        if(t == 1)
            if(exist(out_folder, 'dir'))
                rmdir(out_folder, 's')
            end
            mkdir(out_folder)
        end

        export_fig(fig, out_fig_path, '-nocrop');
    elseif(do_save_fig == 2)
        my_frame = getframe(gcf);
        writeVideo(v, my_frame);
        drawnow;
    else
        drawnow limitrate;
    end
end
if(do_save_fig == 2)
    close(v)
end

export_fig(fig, sprintf('out2/profile_evolution_%d_last_frame.png', pulse_option))

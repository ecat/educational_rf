excitation_kspace_example_number = 1;

rf_pulse_path = './pulses/';

rf_pulse_tag_helper = @(x) 'msinc_highres_tbw_8';

rf_pulse_tag = rf_pulse_tag_helper();
[Nt, dT, ~, ~, pulse_shape_complex, max_B1_G] = ...
    load_designed_pulse(rf_pulse_path, rf_pulse_tag, 0);

keypoints = [];
mag_profile_scale = 1;

gz_1 = ones(Nt, 1); 
if(excitation_kspace_example_number == 1)
    keypoints = ones(Nt, 1);
    keypoints = cat(1, keypoints, 2 * ones(Nt / 2, 1));
    pulse_shape_complex = cat(1, pulse_shape_complex, -zeros(Nt/2, 1));
    gz_1 = cat(1, gz_1, -ones(Nt/2, 1));
    frame_skip = 20; % control frame rate and size of file
elseif(excitation_kspace_example_number == 2)
    % constructed pulse example
    keypoints = ones(Nt, 1);
    keypoints = cat(1, keypoints, 2 * ones(Nt / 2, 1));    
    pulse_shape_complex = cat(1, pulse_shape_complex, -flip(pulse_shape_complex(end/2 + 1:end)));
    gz_1 = cat(1, gz_1, -ones(Nt/2, 1));
    frame_skip = 20;
elseif(excitation_kspace_example_number == 3 || excitation_kspace_example_number == 4)
    % binomial pulse example
    % shave off 1 so never hit index 0 in excitation kspace
    pulse_shape_complex = cat(1, pulse_shape_complex, zeros(Nt - 1, 1), pulse_shape_complex);
    pulse_shape_complex = cat(1, pulse_shape_complex, -zeros(Nt/2, 1));

    if(excitation_kspace_example_number == 4)
        %N_cycles = .75;
        N_cycles = 2.25;
        pulse_shape_complex(1:3 * Nt) = pulse_shape_complex(1: 3 * Nt) .* exp(1i * linspace(0, 2 * pi * N_cycles, 3 * Nt)');
        mag_profile_scale = .5;
    end
    gz_1 = cat(1, gz_1, -gz_1(1:end-1), gz_1);
    gz_1 = cat(1, gz_1, -ones(Nt/2, 1));
    frame_skip = 80;
else
    error();
end

% assumes everything is on-grid (constant gradient with amplitude 1)
Nt_total = numel(gz_1);
excitation_kspace_location_t = cumsum(gz_1);
Nkz = max(excitation_kspace_location_t);
excitation_kspace_value_zt = zeros(Nkz, Nt_total);

for tt = 1:Nt_total
    if(tt > 1)
        excitation_kspace_value_zt(:, tt) = excitation_kspace_value_zt(:, tt - 1);
        excitation_kspace_value_zt(excitation_kspace_location_t(tt), tt) = excitation_kspace_value_zt(excitation_kspace_location_t(tt), tt - 1) + pulse_shape_complex(tt);
    else
        excitation_kspace_value_zt(excitation_kspace_location_t(tt), tt) = pulse_shape_complex(tt);
    end
end

scale = max(abs(excitation_kspace_value_zt(:)));
excitation_kspace_value_zt = excitation_kspace_value_zt / scale;
pulse_shape_complex = pulse_shape_complex / max(abs(pulse_shape_complex));


% need to evaluate excitation kspace at different centers
magnetization_profile_zt = zeros(Nkz, Nt_total);
for tt = 1:Nt_total
     % subtract Nkz/2 so that it is centered when gradient rewinder is done
    excitation_kspace_center = excitation_kspace_location_t(tt) - round(Nkz/2);
    magnetization_profile_zt(:, tt) = cfft(circshift(excitation_kspace_value_zt(:, tt), excitation_kspace_center));
end

%%

crop_height_vals = 1:350;
crop_width_vals = 50:950;
if(1)
    video_writer_helper = @ (v, fig_handle)  writeVideo(v, crop_getframe(getframe(fig_handle), crop_height_vals, crop_width_vals));
else
    video_writer_helper = @(a, b) 1/1;
end

show_complex_channel = excitation_kspace_example_number == 4;

t_to_animate = 1:frame_skip:Nt_total;

fig_excitation_kspace = figure('Position', [100 100 1000 350], 'Color', 'white');
zero_line_color = [.8 .8 .8];
t_undersampled = 1:10:Nt_total;
ts = 1:Nt_total;

subplot(1, 3, 1)
plot([0, Nt], [0 0], 'Color', zero_line_color); hold on;
h1real = animatedline('MaximumNumPoints', Nkz, 'Color', 'k', 'LineWidth', 2);
h1imag = animatedline('MaximumNumPoints', Nkz, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
h1a = animatedline('MaximumNumPoints', 1, 'MarkerSize', 8, 'Marker', 'o', 'MarkerFaceColor', '#8C1515');
h1b = animatedline('MaximumNumPoints', 1, 'MarkerSize', 8, 'Marker', 'o', 'MarkerFaceColor', 'b');
ylim([-1.05 * max(abs(excitation_kspace_value_zt(:))) 1.05 * max(abs(excitation_kspace_value_zt(:)))])


xlim([0, Nt])
xlabel({'Excitation k-Space [1/cm]'}, 'FontSize', 16)
ylabel('Energy [a.u.]', 'FontSize', 16)
xticklabels({})

subplot(1, 3, 2)
mx_line = animatedline('MaximumNumPoints', Nkz, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
my_line = animatedline('MaximumNumPoints', Nkz, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-.');
mxy_line = animatedline('MaximumNumPoints', Nkz, 'Color', 'r', 'LineWidth', 2);
legend('Mx', 'My', '$|$Mxy$|$', 'Location', 'south east', 'FontSize', 10);

ylim([min(real(magnetization_profile_zt(:))) max(abs(magnetization_profile_zt(:)))])
xlim([-Nt/300, Nt/300])
xlabel({'Slice Dimension [cm]'}, 'FontSize', 16)
ylabel('Magnetization [a.u.]', 'FontSize', 16)
xticklabels({})

subplot(2, 3, 3)
plot([0, Nt_total], [0 0], 'Color', zero_line_color); hold on; 
plot(ts(t_undersampled), (real(pulse_shape_complex(t_undersampled))), 'k');
if(show_complex_channel)
    plot(ts(t_undersampled), (imag(pulse_shape_complex(t_undersampled))), 'k:');
end
ylim([-1.1, 1.1])
h2real = animatedline('MaximumNumPoints', 1, 'MarkerSize', 8, 'Marker', 'o', 'MarkerFaceColor', '#8C1515');
h2imag = animatedline('MaximumNumPoints', 1, 'MarkerSize', 8, 'Marker', 'o', 'MarkerFaceColor', 'b');
ylabel('Pulse [a.u.]', 'FontSize', 16)
xticklabels({})

subplot(2, 3, 6)
plot([0, Nt_total], [0 0], 'Color', zero_line_color); hold on; 
plot(ts(t_undersampled), gz_1(t_undersampled), 'k')
ylim([-1.1, 1.1])
h3 = animatedline('MaximumNumPoints', 1, 'MarkerSize', 8, 'Marker', 'o', 'MarkerFaceColor', 'k');
xlabel('Time', 'FontSize', 16)
ylabel('$G_z$ [a.u.]', 'FontSize', 16)
xticklabels({})
%%

keypoint_index = 0;
curr_keypoint_index = 0;
for tt = t_to_animate

    if(numel(keypoints) == 0)
        keypoint_index = 1;
    else
        keypoint_index = keypoints(tt);
    end
    if(keypoint_index > curr_keypoint_index)
        curr_keypoint_index = keypoint_index;
        if(keypoint_index > 1)
            close(v)
        end
        v = VideoWriter(sprintf('out2/excitation_kspace_example_number_%d_step_%d.avi', ...
            excitation_kspace_example_number, keypoint_index), 'Motion JPEG AVI');
        v.Quality = 75;
        v.FrameRate = 30;
        open(v)        
    end

    zlocs = -Nkz/2:Nkz/2-1;
    mx_line.addpoints(zlocs, real(magnetization_profile_zt(:, tt)) * mag_profile_scale)
    my_line.addpoints(zlocs, imag(magnetization_profile_zt(:, tt)) * mag_profile_scale)
    mxy_line.addpoints(zlocs, abs(magnetization_profile_zt(:, tt)) * mag_profile_scale)
    % flip these because excitation kspace reference is flipped
    h1real.addpoints(flip(1:Nkz), real(excitation_kspace_value_zt(:, tt)));
    h1a.addpoints(Nkz - excitation_kspace_location_t(tt), real(excitation_kspace_value_zt(excitation_kspace_location_t(tt), tt)));
    h2real.addpoints(tt, real(pulse_shape_complex(tt)))

    if(show_complex_channel)
        h1imag.addpoints(1:Nkz, imag(excitation_kspace_value_zt(:, tt)));
        h1b.addpoints(excitation_kspace_location_t(tt), imag(excitation_kspace_value_zt(excitation_kspace_location_t(tt), tt)));
        h2imag.addpoints(tt, imag(pulse_shape_complex(tt)))
    end
    h3.addpoints(tt, gz_1(tt))
    drawnow;
    video_writer_helper(v, fig_excitation_kspace);
    %writeVideo(v, getframe(fig_excitation_kspace));
end


close(v)

export_fig(fig_excitation_kspace, sprintf('out2/excitation_kspace_example_number_%d_last_frame.png', excitation_kspace_example_number))

%%

function [ax] = animate_mag_matrix( mag_matrix, rotation_axis_phases, rotation_axis_colorstr, do_drawnow, ...
    time_hist_struct, colors_for_each_spin, varargin)
% mag_matrix, an N_spins x 3 x (t) matrix and plots the spin dimension in one plot
% if the 't' dimension exists, the time course is plotted as a dotted line
% and only the final magnetization is shown
% rotation_axis_phase, an N_spins x N_axesyouwanttodraw vector that is the rotation axis in degrees
%
% time_hist_struct is a structure with 2 variables, can be empty to use gray for the entire timecourse
% Fields are:
% checkpoints, which is the indices of where to change magnetization history
% checkpoint_colors, cell array colour to use, must be same size as checkpoints

if(nargin > 6)
    do_draw_circle = varargin{1};
else
    do_draw_circle = 1;
end

assert(size(mag_matrix, 2) == 3);
assert(size(mag_matrix, 1) == size(rotation_axis_phases, 1) ||numel(rotation_axis_phases) == 0);

do_checkpoint_time_hist = numel(time_hist_struct) ~= 0;
if(do_checkpoint_time_hist)
    assert(numel(time_hist_struct.checkpoints) == numel(time_hist_struct.checkpoint_colors));
end

N_spins = size(mag_matrix, 1);
N_timepoints = size(mag_matrix, 3);

ax = subplot(1, 1, 1); 

% draw coordinate axes    
coord_axes_linewidth = 2;
line([0 0], [0 0], [1 0], 'Color', 'k', 'LineWidth', coord_axes_linewidth);
line([0 0], [1 0], [0 0], 'Color', 'k', 'LineWidth', coord_axes_linewidth);
line([1 0], [0 0], [0 0], 'Color', 'k', 'LineWidth', coord_axes_linewidth);
line([0 0], [0 0], [-1 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', coord_axes_linewidth);
line([0 0], [-1 0], [0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', coord_axes_linewidth);
line([-1 0], [0 0], [0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', coord_axes_linewidth);
set(ax, 'YDir', 'reverse')
text(1.1, 0, 0, 'x', 'HorizontalAlignment', 'left', 'FontSize', 20)
text(0, 1.1, -.05, 'y', 'HorizontalAlignment', 'left', 'FontSize', 20)
text(0, .05, 1.15, 'z', 'HorizontalAlignment', 'left', 'FontSize', 20)

box off;

view(10, 30)
xlim([-1, 1]); 
ylim([-1, 1]); 
zlim([-1, 1]); 

daspect([1 1 1]);
grid on;

hold on;

% plot the xy plane
if(do_draw_circle)
    circle_x = cosd(linspace(0, 360, 360));
    circle_y = sind(linspace(0, 360, 360));
    line(circle_x, circle_y, zeros(360, 1), 'LineStyle', ':', 'LineWidth', .9, 'Color', 'k');
end

for ii = 1:N_spins    
    % plot the magnetization vector
    if( numel(colors_for_each_spin) > 0)
        spin_color = colors_for_each_spin{ii};
    else
        spin_color = 'f3'; %firebrickred
    end
    arrow3([0, 0, 0], mag_matrix(ii, :, end), spin_color, 2, 2);
    
    % plot the rotation axis
    for jj = 1:size(rotation_axis_phases, 2)
        rotation_axis_matrix = .85 * [cosd(rotation_axis_phases(ii, jj)), sind(rotation_axis_phases(ii, jj)), zeros(size(rotation_axis_phases(ii, jj)))];
        arrow3([0, 0, 0], rotation_axis_matrix, rotation_axis_colorstr{jj}, 2, 4);
    end
    
    axis off;
    
    % plot the time history
    if(~do_checkpoint_time_hist)
        hist_x = squeeze(mag_matrix(ii, 1, 1:end-1));
        hist_y = squeeze(mag_matrix(ii, 2, 1:end-1));
        hist_z = squeeze(mag_matrix(ii, 3, 1:end-1));
        line(hist_x, hist_y, hist_z, 'Color', [.3 .3 .3], 'LineStyle', '-', 'LineWidth', .5);
    else
        for t = 1:numel(time_hist_struct.checkpoints)                    
            if(t == 1)
                time_start = 1;
            else
                time_start = time_hist_struct.checkpoints(t - 1);
            end
            time_start = min(time_start, N_timepoints);
            time_end = min(time_hist_struct.checkpoints(t), N_timepoints);
            
            hist_x = squeeze(mag_matrix(ii, 1, time_start:time_end));
            hist_y = squeeze(mag_matrix(ii, 2, time_start:time_end));
            hist_z = squeeze(mag_matrix(ii, 3, time_start:time_end));
            line(hist_x, hist_y, hist_z, 'Color', time_hist_struct.checkpoint_colors{t}, 'LineStyle', ':', 'LineWidth', 2.);
            %line([0, hist_x(1)], [0, hist_y(1)], [0, hist_z(1)], 'Color', [.635, .078, 0.184], 'LineStyle', '-', 'LineWidth', 1)

            arrow3([0, 0, 0], [hist_x(1), hist_y(1), hist_z(1)], 'l', 2, 4);
        end
    end
    
    if(do_drawnow) drawnow; end;
end



end


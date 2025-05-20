function [x, y] = draw_circle_segment(r, theta_min, theta_max)

N = 360;

thetas = linspace(theta_min, theta_max, N)';
is_whole_disc = abs(theta_max - theta_min) > 360;
if is_whole_disc
    x = [];
    y = [];
else
    x = [0];
    y = [0];
end
x = cat(1, x, r * cosd(thetas));
y = cat(1, y, r * sind(thetas));

if ~is_whole_disc
    x(end + 1) = 0;
    y(end + 1) = 0;
end


end


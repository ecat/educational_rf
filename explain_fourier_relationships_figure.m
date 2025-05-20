addPaths()

%multiplication of sinc by cosine
N = 256;
x1 = zeros(N, 1);
x1(end/2 - N/8:end/2 + N/8) = 1;

y2 = cos(linspace(-6 * pi, 6 * pi, N));

y1 = sinc(linspace(-4, 4, N));

x2 = zeros(N, 1);
x2(N/4) = .5;
x2(3 * N / 4) = .5;

x3 = conv(x2, x1, 'same');

title_fontsize = 16;
fig = figure('Position', [100 100 800 400], 'Color', 'white'); 

has = tight_subplot(2, 3, [.2 .1], [.1 .1], [.1 .1]);
axes(has(4))
plot(x1, 'LineWidth', 1.5, 'Color', 'k');
xlim([0, N])
ylim([-1.1, 1.1])
%text(N/2, 1.4, '$\mathrm{rect}(f)$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$*$', 'fontsize', 28, 'HorizontalAlignment', 'center')
text(-1, 1.0, '1', 'FontSize', 20)

axes(has(5))
%stem(x2, 'Color', 'k');
plot(zeros(N, 1), 'Color', 'k'); hold on;
plot([N/4 N/4], [0 x2(N/4)], 'LineWidth', 1.5, 'Color', 'k');
plot(N/4, x2(N/4), 'LineWidth', 1.5, 'Color', 'k', 'Marker', '^');
plot([N*3/4 N*3/4], [0 x2(N/4)], 'LineWidth', 1.5, 'Color', 'k');
plot(N*3/4, x2(N*3/4), 'LineWidth', 1.5, 'Color', 'k', 'Marker', '^');
xlim([0, N])
ylim([-1.1, 1.1])

%text(N/2, 1.4, '$\frac{1}{2}(\delta(f - f_0) + \delta(f + f_0))$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$=$', 'fontsize', 28, 'HorizontalAlignment', 'center')

axes(has(6));
plot(x3, 'LineWidth', 1.5, 'Color', 'k');
xlim([0, N])
ylim([-1.1, 1.1])
text(-5, .5, '$\frac{1}{2}$', 'FontSize', 20)

%text(N/2, 1.4, '$\frac{1}{2}(\mathrm{rect}(f - f_0) + \mathrm{rect}(f + f_0))$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

axes(has(1))
plot(y1, 'k');
xlim([0, N])
ylim([-1.1, 1.1])

text(N/2, 1.4, '$\mathrm{sinc}(t)$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$\odot$', 'fontsize', 28, 'HorizontalAlignment', 'center')

axes(has(2))
plot(real(y2), 'LineWidth', 1.5, 'Color', 'k');
xlim([0, N])
ylim([-1.1, 1.1])
text(N/2, 1.4, '$\cos{(2\pi f_0t)}$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$=$', 'fontsize', 28, 'HorizontalAlignment', 'center')
axes(has(3))
plot(y1 .* y2, 'LineWidth', 1.5, 'Color', 'k'); hold on;
plot(y1, '--', 'LineWidth', 1.5, 'Color', 'r')
xlim([0, N])
ylim([-1.1, 1.1])

text(N/2, 1.4, '$\mathrm{sinc}(t)\cdot \cos{(2\pi f_0t)}$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);
for ii = 1:numel(has)
    axis(has(ii), 'off');
end


   export_fig(fig, sprintf('out2/multiband_explain.png'), '-nocrop');

%%

figb = figure('Position', [100 100 800 400], 'Color', 'white'); 

x1 = zeros(N, 1);
x1(end/2 - N/32:end/2 + N/32) = 1;
x2 = zeros(N, 1);
x2_1_indices = 16:32:N;
x2(x2_1_indices) = 1;

y2 = zeros(N, 1);
y2_1_indices = 32:64:N;
y2(y2_1_indices) = 1;

x3 = conv(x1, y2, 'same');

has = tight_subplot(2, 3, [.2 .1], [.1 .1], [.1 .1]);
axes(has(4))
plot(x1, 'LineWidth', 1.5, 'Color', 'k');
xlim([0, N])
ylim([-1.1, 1.1])
%text(N/2, 1.4, '$\mathrm{rect}(f)$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$*$', 'fontsize', 28, 'HorizontalAlignment', 'center')


axes(has(5))
%stem(x2, 'Color', 'k');
plot(zeros(N, 1), 'Color', 'k'); hold on;
for xx = y2_1_indices
    plot([xx xx], [0 y2(xx)], 'LineWidth', 1.5, 'Color', 'k');
    plot(xx, y2(xx), 'LineWidth', 1.5, 'Color', 'k', 'Marker', '^');
end
xlim([0, N])
ylim([-1.1, 1.1])

%text(N/2, 1.4, '$\frac{1}{2}(\delta(f - f_0) + \delta(f + f_0))$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$=$', 'fontsize', 28, 'HorizontalAlignment', 'center')

axes(has(6));
plot(x3, 'LineWidth', 1.5, 'Color', 'k');
xlim([0, N])
ylim([-1.1, 1.1])
%text(-5, .5, '$\frac{1}{2}$', 'FontSize', 20)

%text(N/2, 1.4, '$\frac{1}{2}(\mathrm{rect}(f - f_0) + \mathrm{rect}(f + f_0))$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

axes(has(1))
plot(y1, 'k');
xlim([0, N])
ylim([-1.1, 1.1])

text(N/2, 1.4, '$\mathrm{sinc}(t)$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$\odot$', 'fontsize', 28, 'HorizontalAlignment', 'center')

axes(has(2))
plot(zeros(N, 1), 'Color', 'k'); hold on;
for xx = x2_1_indices
    plot([xx xx], [0 x2(xx)], 'LineWidth', 1.5, 'Color', 'k'); hold on;
    plot(xx, x2(xx), 'LineWidth', 1.5, 'Color', 'k', 'Marker', '^');
end
xlim([0, N])
ylim([-1.1, 1.1])
text(N/2, 1.4, '$\Sigma_k{\delta(t - kT_0)}$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize);

text(1.25 * N, 0, '$=$', 'fontsize', 28, 'HorizontalAlignment', 'center')
axes(has(3))
plot(y1 .* x2.', 'LineWidth', 1.5, 'Color', 'k'); hold on;
plot(y1, '--', 'LineWidth', 1.5, 'Color', 'r')
xlim([0, N])
ylim([-1.1, 1.1])

text(N/2, 1.4, '$\mathrm{sinc}(t)\cdot \Sigma_k{\delta(t - kT_0)}$', 'HorizontalAlignment', 'center', 'FontSize', title_fontsize, 'Interpreter', 'latex');
for ii = 1:numel(has)
    axis(has(ii), 'off');
end


export_fig(figb, sprintf('out2/discrete_sampling_explain.png'), '-nocrop');

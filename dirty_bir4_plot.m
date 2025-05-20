[mag, phi] = get_bir4_pulse(.18, 1000, .05, 2, 90, 10, 1);
ts = linspace(0, 6, 1000);
figure('Position', [100 100 350 600], 'Color', 'white');
subplot(211);
plot(ts, mag * 100, 'LineWidth', 2);
ylim([0 22])
ylabel('Magnitude [uT]')
subplot(212);
plot(ts, mod(phi, 2 * pi) - pi, 'LineWidth', 2)
ylabel('Phase [rad]')
ylim([-4 4])


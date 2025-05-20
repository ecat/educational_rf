function [N, dT, dfs, M, pulse_shape_complex, max_B1_G] = load_designed_pulse(rf_pulse_path, rf_pulse_tag, do_plot_pulse)


suffixes = {'pulse', 'dT', 'dfs', 'M'};

get_filepath = @(tag, suffix) strcat(rf_pulse_path, rf_pulse_tag, '_', suffix, '.txt');
loaded_files = {};

for suffix = suffixes

    path_to_load = get_filepath(rf_pulse_tag, suffix{1});

    if(exist(path_to_load, 'file'))
        loaded_files{end+1} = load(path_to_load);
    end
end

if(size(loaded_files{1}, 2) == 2)
    pulse_shape_G = loaded_files{1};
else
    pulse_shape_G = loaded_files{1}.';
end
assert(size(pulse_shape_G, 2) == 2, 'pulse shape should be saved as real channel and imag channel');
N = size(pulse_shape_G, 1);
dT = loaded_files{2};
duration = N * dT;

has_simulation = numel(loaded_files) > 2;

if(has_simulation)
    dfs = loaded_files{3};
    M = loaded_files{4};
else
    dfs = 0;
    M = 0;
end


pulse_shape_complex = pulse_shape_G(:, 1) + 1i * pulse_shape_G(:, 2);
max_B1_G = max(abs(pulse_shape_complex));

if do_plot_pulse 
    pulse_shape_G = cat(2, real(pulse_shape_complex), imag(pulse_shape_complex));
    figure; plot(linspace(0, duration, N), (pulse_shape_G).'); legend('re', 'im'); ylabel('G'); xlabel('T [ms]');
end

end


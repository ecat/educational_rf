function [M_result, G] = simulate_rfpulse_every_timepoint(dfs, B1_t, ...
    N_timepoints, dT, T_rewinder, N_rewinder_timepoints, G_in, constant_off_resonance, varargin)

N_spins = numel(dfs);


if(nargin > 8)
    M_init_z = varargin{1};
else
    M_init_z = repmat([0; 0; 1], [1 N_spins]);
end

if(numel(G_in) == 0)
    M_result = zeros(3, N_spins, N_timepoints + N_rewinder_timepoints);
    G = ones(size(M_result, 3), 1);
    G(N_timepoints + 1:end) = -T_rewinder / (N_rewinder_timepoints * dT);
else
    assert(numel(G_in) == numel(B1_t));
    assert(T_rewinder == 0, 'if G is provided assume includes rewinder')
    assert(N_rewinder_timepoints == 0)
    %assert(all(abs(G_in) - 1 < eps))

    M_result = zeros(3, N_spins, N_timepoints);
    G = G_in;
end


gamma = 4.258; % kHz / G

for ss = 1:N_spins
    M_next = M_init_z(:, ss);
    for tt = 1:N_timepoints
        M_next = throt(abs(B1_t(tt)) * dT * gamma * 360, angle(B1_t(tt)) * 180 / pi) * M_next;

        if(numel(G_in) == 0)
            M_next = zrot((dfs(ss) + constant_off_resonance) * dT * 360) * M_next;
        else
            M_next = zrot((G_in(tt) * dfs(ss) + constant_off_resonance)* dT * 360) * M_next;
        end
        
        M_result(:, ss, tt) = M_next;
    end
    
    for tt = 1:N_rewinder_timepoints
        M_next = zrot((-dfs(ss) + constant_off_resonance) * T_rewinder / N_rewinder_timepoints * 360) * M_next;
        M_result(:, ss, tt + N_timepoints) = M_next;
    end
end


end

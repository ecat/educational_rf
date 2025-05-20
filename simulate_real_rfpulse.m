function [M_result_rzt] = simulate_real_rfpulse(dfs, B1_t, dT) 
% dfs [kHz]
% B1_t [G]
% dT [ms]
% N_timepoints integer
% returns magnetization profile for each timepoint of B1_t
assert(isvector(B1_t), 'B1_t should be 1D') 
assert(isreal(B1_t), 'This function assumes B1_t is real')

N_timepoints = numel(B1_t);
N_spins = numel(dfs);
M_result_rzt = zeros(3, N_spins, N_timepoints);

gamma = 4.258; % [kHz / G]

for rr = 1:N_spins % loop over locations in slice
    M_next = [0; 0; 1]; % initialize at equilibrium
    for tt = 1:N_timepoints
        Ry = yrot(B1_t(tt) * gamma * dT * 360); % RF on, gradient off
        Rz = zrot(dfs(rr) * dT * 360); % RF off, gradient on
        M_next = Rz * Ry * M_next;
        M_result_rzt(:, rr, tt) = M_next;
    end
end

% http://hilite.me/ perldoc

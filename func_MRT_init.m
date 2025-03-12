function [Theta,G_matrices,debug] = func_MRT_init(H_RI,H_IT,properties,Sl_blocks)
c = 3e8;
lambda = c/properties.carrier_freq;
N = properties.No_SIM_RE_z*properties.No_SIM_RE_y;
L = properties.No_SIM_layers;
K = properties.No_Users;

Sl_21 = properties.Sl_blocks.S21;
[Theta{L},obj_f] = func_MRT_GC(H_IT,H_RI,1);
Tl = properties.Tl;

H_temp = H_RI*Theta{L};
for i = L-1 :-1: 1
    H_temp = H_temp*Sl_21;
    [Theta{i},obj_f] = func_MRT_GC(H_IT,H_temp,1);
    H_temp = H_temp*Theta{i};
end
debug.H_wrong = H_temp*H_IT;
for l = 1:L
    G_matrices{l} = blkdiag(Theta{l}, conj(Theta{l}));
end

T_I = G_matrices{1} * Tl;
for l = 2:L-1
    T_I = T_I * G_matrices{l} * Tl;
end
T_I = T_I * G_matrices{L};
% Extract T_I,22
T_I_22 = T_I(N+1:end, N+1:end);

% Compute H
debug.H_correct = H_RI * inv(T_I_22) * H_IT;
debug.H_correct_loss = norm(debug.H_correct - diag(diag(debug.H_correct)),'fro')^2;


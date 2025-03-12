function [Theta,R] = func_MRT_GC(W_norm,H_norm,NG)
% Compute the scattering matrix Theta given the normalized channels
% W_norm and H_norm, and the group size NG
% Inputs:   H_norm: normalized channel H
%           W_norm: normalized channel W
%           NG: group size (NG=N means fully connected)
%                          (NG=1 means single connected)
% Outputs:  Theta: scattering matrix symuni
%           R: scattering matrix relaxed 
N = size(W_norm,1);
G_matrix = W_norm*H_norm;
Theta = zeros(N,N);
R = zeros(N,N);
if(NG == 1) % Single connected
    R = diag(G_matrix');
    Theta = diag(R./abs(R));
    R = diag(R);
else % Group or fully connected
end
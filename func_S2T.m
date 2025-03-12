% This function transforms S-parameters to T-parameters
% S has to be a square of size 2N x 2N
function [T,T_TT,T_TR,T_RT,T_RR] = func_S2T(S)
[n m] = size(S);
if (n==m)
    if (mod(n,2)==0)
    N = n/2;
    else
        msg = 'Size of S must be even.';
        error(msg)
    end
else
    msg = 'S has to be square.';
    error(msg)
end

S_TT = S(1:N,1:N);
S_TR = S(1:N,N+1:2*N);
S_RT = S(N+1:2*N,1:N);
S_RR = S(N+1:2*N,N+1:2*N);

T_TT = S_TR - S_TT*inv(S_RT)*S_RR;
T_TR = S_TT*inv(S_RT);
T_RT = -inv(S_RT)*S_RR;
T_RR = inv(S_RT);
T = [T_TT T_TR;...
     T_RT T_RR];
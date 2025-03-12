function [Sl,Zl,Sl_blocks,Zl_blocks,debug] = func_SIM_MC_dipole(properties)
% clc
% clear
% close all
% properties.carrier_freq = 28e9;
% properties.No_SIM_RE_z = 2;
% properties.No_SIM_RE_y = 2;
% properties.No_SIM_layers = 3;
c = 3e8;                           % Speed of light
% epsilon0 = 8.8541878188e-12;       % Vacuum electric permittivity
% mu0 = 1.25663706127e-6;            % Vacuum magnetic permeability
% eta0 = sqrt(mu0/epsilon0);         % Vacuum characteristic impedance
eta0 = 377;
f = properties.carrier_freq;       % Carrier frequency
lambda = c/f;                      % Wavelength
k0 = 2*pi/lambda;                  % Wave number
% properties.RE_length = lambda/4;
% properties.RE_raduis = lambda/500;

L = properties.No_SIM_layers;   % No. SIM layers
Ny = properties.No_SIM_RE_y;    % RE in y-axis
Nz = properties.No_SIM_RE_z;    % RE in z-axis
N = Ny*Nz;                      % RE in y-z plane
T_SIM = (L-1)*properties.Inter_Layer_Sep;
L_xtot = T_SIM;                       % SIM Thickness
l_x = L_xtot/(L-1);                   % Adjacent SIM layers seperation
l_z = properties.Adjac_Elem_Sep_z;    % Adjacent REs seperation in z
l_y = properties.Adjac_Elem_Sep_y;    % Adjacent REs seperation in y
L_z = properties.RE_length;           % RE size in z
L_y = properties.RE_raduis;           % RE size in y
L_ztot = L_z*Nz + (Nz - 1)*l_z;       % SIM size in z
L_ytot = L_y*Ny + (Ny - 1)*l_y;       % SIM size in y

% lth Layer REs coordinates
for i = 1:L
    for m = 1:Ny
        for n = 1:Nz
            x_sn{i}(m,n) = l_x*(i-1);
            y_sn{i}(m,n) = (m - 1).*l_y;
            z_sn{i}(m,n) = (n - 1).*l_z;
            r_sn{i}{m,n} = [x_sn y_sn z_sn];
        end
    end
    x_sn{i} = x_sn{i}(:);
    y_sn{i} = y_sn{i}(:);
    z_sn{i} = z_sn{i}(:);
end

% Mutual coupling between the RIS elements
Zl11 = zeros(N,N);
for a = 1:N
    for b = 1:N
        if a == b
            d_xy11(a,b) = L_y^2;
        else
            d_xy11(a,b) = (x_sn{1}(a) - x_sn{1}(b)).^2 + (y_sn{1}(a) - y_sn{1}(b)).^2;
        end
        fun = @(z1,z2) cost(z1,z2,d_xy11(a,b),z_sn{1}(a),z_sn{1}(b));
        Zl11(a,b) = integral2(fun, z_sn{1}(a)-L_z/2, z_sn{1}(a)+L_z/2,...
                    z_sn{1}(b)-L_z/2, z_sn{1}(b)+L_z/2);
    end
end

Zl22 = Zl11;

Zl12 = zeros(N,N);
for a = 1:N
    for b = 1:N
        d_xy12(a,b) = (x_sn{1}(a) - x_sn{2}(b)).^2 + (y_sn{1}(a) - y_sn{2}(b)).^2;
        fun = @(z1,z2) cost(z1,z2,d_xy12(a,b),z_sn{1}(a),z_sn{2}(b));
        Zl12(a,b) = integral2(fun, z_sn{1}(a)-L_z/2, z_sn{1}(a)+L_z/2,...
                    z_sn{2}(b)-L_z/2, z_sn{2}(b)+L_z/2);
    end
end

Zl21 = Zl12.';

Z0 = 50;
Zl = [Zl11 Zl12;...
      Zl21 Zl22];
Zl(~~eye(2*N,2*N)) = 50;
Sl = (Zl + Z0*eye(2*N,2*N))\(Zl - Z0*eye(2*N,2*N));

Sl_blocks.S11 = Sl(1:N,1:N);
Sl_blocks.S12 = Sl(1:N,N+1:2*N);
Sl_blocks.S21 = Sl(N+1:2*N,1:N);
Sl_blocks.S22 = Sl(N+1:2*N,N+1:2*N);

% Sl_blocks.Sl11 = (Zl11 + Z0)\(Zl11 - Z0);
% Sl_blocks.Sl12 = (Zl12 + Z0)\(Zl12 - Z0);
% Sl_blocks.Sl21 = (Zl21 + Z0)\(Zl21 - Z0);
% Sl_blocks.Sl22 = (Zl22 + Z0)\(Zl22 - Z0);

Zl_blocks.Z11 = Zl11;
Zl_blocks.Z12 = Zl12;
Zl_blocks.Z21 = Zl21;
Zl_blocks.Z22 = Zl22;
debug.d_xy11 = d_xy11;
debug.d_xy12 = d_xy12;
debug.x_sn = x_sn;
debug.y_sn = y_sn;
debug.z_sn = z_sn;
    function obj = cost(z1,z2,d_xy,z_a,z_b)
        d_ab = sqrt(d_xy + (z1-z2).^2);
        obj = 1j*eta0/(4*pi*k0)...
            .*((z1-z2).^2./d_ab.^2.*(3./d_ab.^2 + 1j*3*k0./d_ab - k0^2) - (1j*k0 + 1./d_ab)./d_ab + k0.^2)...
            .*exp(-1j*k0.*d_ab)./d_ab ...
            .*sin(k0*(L_z/2 - abs(z1 - z_a)))./sin(k0*L_z/2)...
            .*sin(k0*(L_z/2 - abs(z2 - z_b)))./sin(k0*L_z/2);
    end
end

function [Sl_21,Sl_21_norm] = func_SIM_RaySom(properties)

c = 3e8;                           % Speed of light
f = properties.carrier_freq;       % Carrier frequency
lambda = c/f;                      % Wavelength

L = properties.No_SIM_layers;   % No. SIM layers
Ny = properties.No_SIM_RE_y;    % RE in y-axis
Nz = properties.No_SIM_RE_z;    % RE in z-axis
N = Ny*Nz;                      % RE in y-z plane
T_SIM = (L-1)*properties.Inter_Layer_Sep;
L_xtot = T_SIM;                       % SIM Thickness
l_x = L_xtot/(L-1);                   % Adjacent SIM layers seperation
l_z = properties.Adjac_Elem_Sep_z;      % Adjacent REs seperation in z
l_y = properties.Adjac_Elem_Sep_y;      % Adjacent REs seperation in y
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


for a = 1:N
    for b = 1:N
        d_n(a,b) = sqrt((x_sn{1}(a) - x_sn{2}(b)).^2 + (y_sn{1}(a) - y_sn{2}(b)).^2 + ...
            (z_sn{1}(a) - z_sn{2}(b)).^2); % distance from lth SIM layer to (l-1)th SIM layer
    end
end


cos_chi_n = l_x./d_n;

Sl_21 = L_y*L_z.*cos_chi_n./d_n.*(1/2/pi./d_n-1j*1/lambda).*exp(1j*2*pi*d_n./lambda);
Sl_21_norm = norm(Sl_21,'fro').^2;

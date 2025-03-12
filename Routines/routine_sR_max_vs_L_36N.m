clc
clear
close all
properties.iso = false;
c = 3e8; 
properties.carrier_freq = 28e9;
lambda = c/properties.carrier_freq;                      
T_SIM = lambda/2/6;
properties.RE_length = lambda/4;
properties.RE_raduis = lambda/500;
properties.No_Users = 2;
K = properties.No_Users;
M = K;
properties.GDA.maxIter = 600;
properties.GDA.tol = 0.005/100;
properties.GDA.alpha = 1; % Step size
properties.GDA.delta = 1e-6;
properties.GDA.tau = 0.8;

Iter = 100;

%%
L_range = [2:4,6];
N_range = ones(1:4).*(36);
H_RI_all = sqrt(1)./sqrt(2)*(randn(K,N_range(1),Iter) + 1j*randn(K,N_range(1),Iter));
H_IT_all = sqrt(1)./sqrt(2)*(randn(N_range(1),K,Iter) + 1j*randn(N_range(1),K,Iter));
%%
P_tx = 1;
N0 = 1;
for i_loop = 1:length(L_range)
    i_loop
    tic
    warning('off','all')
    iter = 0;
    r_ee = 0;
    r_se = 0;
    r_ss = 0;
    L = L_range(i_loop);
    N = N_range(i_loop);
    properties.No_SIM_layers = L;
    properties.Inter_Layer_Sep = (T_SIM)./(L-1);
    properties.No_SIM_RE_z = 6;
    properties.No_SIM_RE_y = N./properties.No_SIM_RE_z;
    properties.Adjac_Elem_Sep_z = lambda/2;
    properties.Adjac_Elem_Sep_y = lambda/2./(N./36);
    [Sl,Zl,Sl_blocks,Zl_blocks] = func_SIM_MC_dipole(properties);
    Sl_21 = Sl_blocks.S21;
    [Tl] = func_S2T(Sl);
    properties.Sl_blocks.S21 = Sl_21;
    properties.Tl = Tl;
    P_max = P_tx;
    while(iter<Iter)
        iter
        H_IT = H_IT_all(1:N,1:M,iter+1);
        H_RI = H_RI_all(1:K,1:N,iter+1);
        [~,opt_ee] = func_sR_MAX_GDA_ExactExact(H_RI,H_IT,Tl,N0,properties);
        [~,opt_se] = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);
        r_ee = r_ee + func_compute_sR(eye(K,K),opt_ee.H_opt,N0);
        r_se = r_se + func_compute_sR(eye(K,K),opt_se.H_opt,N0);
        r_ss = r_ss + func_compute_sR(eye(K,K),opt_se.H_opt_ss,N0);
        iter = iter + 1;
    end
    R_ee(i_loop) = r_ee./Iter;
    R_se(i_loop) = r_se./Iter;
    R_ss(i_loop) = r_ss./Iter;
end
%%
figure
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 2.5, 3.5], 'PaperUnits', 'Inches', 'PaperSize', [2.5, 3.5]);
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(gcf,'color','w');
t = tiledlayout(1,1,'TileSpacing','compact');
% t.Padding = 'compact';
ax1 = axes(t);
ax1.ColorOrderIndex = 1;
% plot(ax1,L_range,R_ee,'-s','LineWidth',1.25); hold on
% plot(ax1,L_range,R_se,'-^','LineWidth',1.25); hold on
% plot(ax1,L_range,R_ss,'-*','LineWidth',1.25); hold on
dataMatrix = [R_ee(:), R_se(:), R_ss(:)];
% Create a grouped bar chart. 
bar(ax1, 1:4, dataMatrix, 'grouped');

ax1.XGrid = 'on';
ax1.YGrid = 'on';
xlabel(ax1,'Number of layers','interpreter','tex','fontsize',12)
ylabel(ax1,'Sum rate (bps/Hz)','interpreter','tex','fontsize',12)
legend(ax1, {'EE', 'SE', 'SS'}, 'interpreter', 'tex', 'fontsize', 8, 'location', 'northeast')
ax1.XLim = [1-0.5, 4+0.5];

% legend(ax1,'EE','SE','SS',...
%            'interpreter','tex','fontsize',8,'location','northeast') 
% ax1.YLim = [0 10.5];
ax1.FontSize = 10;
ax1.LineWidth = 0.75;
ax1.XTick = 1:4;
ax1.XTickLabel = {'2', '3', '4', '6'};
grid on
grid minor 
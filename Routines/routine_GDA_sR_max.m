clc
clear
close all
c = 3e8;    
properties.carrier_freq = 28e9;
properties.No_SIM_RE_z = 6;
properties.No_SIM_RE_y = 6;
properties.No_SIM_layers = 3;
lambda = c/properties.carrier_freq;                      
% properties.RE_length = lambda/4;
% properties.RE_raduis = lambda/500;
properties.RE_length = lambda/4;
properties.RE_raduis = lambda/4;

N = properties.No_SIM_RE_z*properties.No_SIM_RE_y;
L = properties.No_SIM_layers;
properties.No_Users = 2;
K = properties.No_Users;
properties.GDA.maxIter = 600;
properties.GDA.tol = 0.005/100;
properties.GDA.alpha = 1; % Step size
properties.GDA.delta = 1e-6;
properties.GDA.tau = 0.8;

N0 = 1;
load routine_GDA_sR_max_saved_channels
%%
properties.Inter_Layer_Sep = lambda/2;
properties.Adjac_Elem_Sep = lambda/2;
properties.Adjac_Elem_Sep_z = properties.Adjac_Elem_Sep;
properties.Adjac_Elem_Sep_y = properties.Adjac_Elem_Sep;
[Sl,Zl,Sl_blocks,Zl_blocks] = func_SIM_MC_dipole(properties);
Sl_21 = Sl_blocks.S21;
[Tl] = func_S2T(Sl);
properties.Sl_blocks.S21 = Sl_blocks.S21;
properties.Tl = Tl;
[sR_A_EE,opt_A_EE] = func_sR_MAX_GDA_ExactExact(H_RI,H_IT,Tl,N0,properties);
[sR_A_SS,opt_A_SE] = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);
[Sl_21,~] = func_SIM_RaySom(properties);
properties.Sl_blocks.S21 = Sl_21;
Sl = zeros(2*N,2*N);
Sl(N+1:end,1:N) = Sl_21;
[Tl] = func_S2T(Sl);
properties.Tl = Tl;
sR_A_SS_RaySom = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);

properties.Inter_Layer_Sep = lambda/3;
properties.Adjac_Elem_Sep = lambda/3;
properties.Adjac_Elem_Sep_z = properties.Adjac_Elem_Sep;
properties.Adjac_Elem_Sep_y = properties.Adjac_Elem_Sep;

[Sl,Zl,Sl_blocks,Zl_blocks] = func_SIM_MC_dipole(properties);
Sl_21 = Sl_blocks.S21;
[Tl] = func_S2T(Sl);
properties.Sl_blocks.S21 = Sl_blocks.S21;
properties.Tl = Tl;
[sR_B_EE,opt_B_EE] = func_sR_MAX_GDA_ExactExact(H_RI,H_IT,Tl,N0,properties);
[sR_B_SS,opt_B_SE] = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);
[Sl_21,Sl_21_norm] = func_SIM_RaySom(properties);
properties.Sl_blocks.S21 = Sl_21;
Sl = zeros(2*N,2*N);
Sl(N+1:end,1:N) = Sl_21;
[Tl] = func_S2T(Sl);
properties.Tl = Tl;
sR_B_SS_RaySom = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);

properties.Inter_Layer_Sep = lambda/3;
properties.Adjac_Elem_Sep = lambda/2;
properties.Adjac_Elem_Sep_z = properties.Adjac_Elem_Sep;
properties.Adjac_Elem_Sep_y = properties.Adjac_Elem_Sep;

[Sl,Zl,Sl_blocks,Zl_blocks] = func_SIM_MC_dipole(properties);
Sl_21 = Sl_blocks.S21;
[Tl] = func_S2T(Sl);
properties.Sl_blocks.S21 = Sl_blocks.S21;
properties.Tl = Tl;
[sR_C_EE,opt_C_EE] = func_sR_MAX_GDA_ExactExact(H_RI,H_IT,Tl,N0,properties);
[sR_C_SS,opt_C_SE] = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);
[Sl_21,Sl_21_norm] = func_SIM_RaySom(properties);
properties.Sl_blocks.S21 = Sl_21;
Sl = zeros(2*N,2*N);
Sl(N+1:end,1:N) = Sl_21;
[Tl] = func_S2T(Sl);
properties.Tl = Tl;
sR_C_SS_RaySom = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);

properties.Inter_Layer_Sep = lambda/2;
properties.Adjac_Elem_Sep = lambda/3;
properties.Adjac_Elem_Sep_z = properties.Adjac_Elem_Sep;
properties.Adjac_Elem_Sep_y = properties.Adjac_Elem_Sep;

[Sl,Zl,Sl_blocks,Zl_blocks] = func_SIM_MC_dipole(properties);
Sl_21 = Sl_blocks.S21;
[Tl] = func_S2T(Sl);
properties.Tl = Tl;
properties.Sl_blocks.S21 = Sl_blocks.S21;
[sR_D_EE,opt_D_EE] = func_sR_MAX_GDA_ExactExact(H_RI,H_IT,Tl,N0,properties);
[sR_D_SS,opt_D_SE] = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);
[Sl_21,Sl_21_norm] = func_SIM_RaySom(properties);
properties.Sl_blocks.S21 = Sl_21;
Sl = zeros(2*N,2*N);
Sl(N+1:end,1:N) = Sl_21;
[Tl] = func_S2T(Sl);
properties.Tl = Tl;
sR_D_SS_RaySom = func_sR_MAX_GDA_SimplExact(H_RI,H_IT,Tl,N0,properties);

%%
figure
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 3.5], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 3.5]);
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(gcf,'color','w');
t = tiledlayout(1,1,'TileSpacing','compact');
% t.Padding = 'compact';
ax1 = axes(t);
ax1.ColorOrderIndex = 1;
plot(ax1,NaN,NaN,'-k','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'--k','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,':k','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,NaN,NaN,'s','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'^','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'o','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'x','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(sR_A_EE)],sR_A_EE(1:1:end),'-s','MarkerIndices',[1:100:length(sR_A_EE)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_B_EE)],sR_B_EE(1:1:end),'-^','MarkerIndices',[1:100:length(sR_B_EE)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_C_EE)],sR_C_EE(1:1:end),'-o','MarkerIndices',[1:100:length(sR_C_EE)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_D_EE)],sR_D_EE(1:1:end),'-x','MarkerIndices',[1:100:length(sR_D_EE)],'LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(opt_A_SE.sum_rate_simpleExact)],opt_A_SE.sum_rate_simpleExact(1:1:end),'--s','MarkerIndices',[1:100:length(opt_A_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(opt_B_SE.sum_rate_simpleExact)],opt_B_SE.sum_rate_simpleExact(1:1:end),'--^','MarkerIndices',[1:100:length(opt_B_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(opt_C_SE.sum_rate_simpleExact)],opt_C_SE.sum_rate_simpleExact(1:1:end),'--o','MarkerIndices',[1:100:length(opt_C_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(opt_D_SE.sum_rate_simpleExact)],opt_D_SE.sum_rate_simpleExact(1:1:end),'--x','MarkerIndices',[1:100:length(opt_D_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(sR_A_SS)],sR_A_SS(1:1:end),':s','MarkerIndices',[1:100:length(sR_A_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_B_SS)],sR_B_SS(1:1:end),':^','MarkerIndices',[1:100:length(sR_B_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_C_SS)],sR_C_SS(1:1:end),':o','MarkerIndices',[1:100:length(sR_C_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_D_SS)],sR_D_SS(1:1:end),':x','MarkerIndices',[1:100:length(sR_D_SS)],'LineWidth',1.25); hold on
ax1.XGrid = 'on';
ax1.YGrid = 'on';
xlabel(ax1,'Number of iterations','interpreter','tex','fontsize',12)
ylabel(ax1,'Sum rate (bps/Hz)','interpreter','tex','fontsize',12)
legend(ax1,'EE','SE','SS',...
           'Case 1','Case 2','Case 3','Case 4',...
           'interpreter','tex','fontsize',8,'location','northeast') 
ax1.XLim = [1 625];
% ax1.YLim = [5e0 1e14];
ax1.FontSize = 10;
ax1.LineWidth = 0.75;

%%
figure
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 3.5], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 3.5]);
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(gcf,'color','w');
t = tiledlayout(1,1,'TileSpacing','compact');
% t.Padding = 'compact';
ax1 = axes(t);
ax1.ColorOrderIndex = 1;
plot(ax1,NaN,NaN,'-k','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'--k','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,':k','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'-.k','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,NaN,NaN,'s','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'^','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'o','LineWidth',1.25); hold on
plot(ax1,NaN,NaN,'x','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(sR_A_EE)],sR_A_EE(1:1:end),'-s','MarkerIndices',[1:100:length(sR_A_EE)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_B_EE)],sR_B_EE(1:1:end),'-^','MarkerIndices',[1:100:length(sR_B_EE)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_C_EE)],sR_C_EE(1:1:end),'-o','MarkerIndices',[1:100:length(sR_C_EE)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_D_EE)],sR_D_EE(1:1:end),'-x','MarkerIndices',[1:100:length(sR_D_EE)],'LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(opt_A_SE.sum_rate_simpleExact)],opt_A_SE.sum_rate_simpleExact(1:1:end),'--s','MarkerIndices',[1:100:length(opt_A_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(opt_B_SE.sum_rate_simpleExact)],opt_B_SE.sum_rate_simpleExact(1:1:end),'--^','MarkerIndices',[1:100:length(opt_B_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(opt_C_SE.sum_rate_simpleExact)],opt_C_SE.sum_rate_simpleExact(1:1:end),'--o','MarkerIndices',[1:100:length(opt_C_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(opt_D_SE.sum_rate_simpleExact)],opt_D_SE.sum_rate_simpleExact(1:1:end),'--x','MarkerIndices',[1:100:length(opt_D_SE.sum_rate_simpleExact)],'LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(sR_A_SS)],sR_A_SS(1:1:end),':s','MarkerIndices',[1:100:length(sR_A_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_B_SS)],sR_B_SS(1:1:end),':^','MarkerIndices',[1:100:length(sR_B_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_C_SS)],sR_C_SS(1:1:end),':o','MarkerIndices',[1:100:length(sR_C_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_D_SS)],sR_D_SS(1:1:end),':x','MarkerIndices',[1:100:length(sR_D_SS)],'LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
plot(ax1,[1:1:length(sR_A_SS_RaySom)],sR_A_SS_RaySom(1:1:end),'-.s','MarkerIndices',[1:100:length(sR_A_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_B_SS_RaySom)],sR_B_SS_RaySom(1:1:end),'-.^','MarkerIndices',[1:100:length(sR_B_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_C_SS_RaySom)],sR_C_SS_RaySom(1:1:end),'-.o','MarkerIndices',[1:100:length(sR_C_SS)],'LineWidth',1.25); hold on
plot(ax1,[1:1:length(sR_D_SS_RaySom)],sR_D_SS_RaySom(1:1:end),'-.x','MarkerIndices',[1:100:length(sR_D_SS)],'LineWidth',1.25); hold on

ax1.XGrid = 'on';
ax1.YGrid = 'on';
xlabel(ax1,'Number of iterations','interpreter','tex','fontsize',12)
ylabel(ax1,'Sum rate (bps/Hz)','interpreter','tex','fontsize',12)
legend(ax1,'EE','SE','SS','SS-RS',...
           'Case 1','Case 2','Case 3','Case 4',...
           'interpreter','tex','fontsize',8,'location','northeast') 
ax1.XLim = [1 625];
% ax1.YLim = [5e0 1e14];
ax1.FontSize = 10;
ax1.LineWidth = 0.75;

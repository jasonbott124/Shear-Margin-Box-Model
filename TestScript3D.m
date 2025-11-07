% Takahe_3D_Heat_Diffusion

%% Intitialize Paramaters
Parameters = build_physical_parameters("enthalpy");
% Surf_ice_Enthalpy = @(Surf_Temp,Parameters) Surf_Temp*Parameters.rho.ice*Parameters.cp.ice;

% %% Define Thermo-Physical Properties
[phi,kappa_sys,rhocp_sys,H,rho,cp,kappa,h,chi] = system_properties_grl(Parameters);
% plot_system_properties(phi,Ts,Tl,rhocp_sys,kappa_sys,col)

%% Define Grid(s) (Discrete Space(s))

% Heat Grid
Grid.xmin = 0; Grid.xmax = 5;
Grid.ymin = 0; Grid.ymax = 1;
Grid.zmin = 0; Grid.zmax = 0.3; 
Grid.Nz = 50; Grid.Ny = 165; Grid.Nx = 165*5;
Grid.geom = 'cartesian';

Grid = build_grid_3D(Grid);
[D,G,~,I] = build_ops_3D(Grid);

%% Boundary Conditions

% Dirchelet BC's
BC.dof_dir = [Grid.dof_zmax];
BC.dof_f_dir = [];
BC.g = -20 * ones(length(Grid.dof_zmax),1); % Surface temp assigned here

% Neumanm BC's
BC.dof_neu   = [Grid.dof_zmin];
BC.dof_f_neu = [Grid.dof_f_zmin];
BC.qb        = 90/1000; % 90 mW/m^2

[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve Heat Diffusion
Kd = 2.2;
% Kd = @(u) comp_mean(reshape(kappa_sys(u),Grid.Ny,Grid.Nx),-1,1,Grid); 
% fs = sparse(Grid.N,1,0);
L = -D*Kd*G;

% Define the full kappa_N(u) as a function
% kappa_N = @(u) set_kappa_field(u, Grid, Geom, kappa, kappa_sys);
% 
% H = @(u) set_h_field(u, Grid, Geom, H, h_sys);
% 
% % Then define Kd using kappa_N
% Kd = @(u) comp_mean(reshape(kappa_N(u), Grid.Ny, Grid.Nx), -1, 1, Grid);
% 
% [uss,Lss,rss,u_h,Kd_h] = compute_geotherm(Kd,D,G,fs,B,N,fn,BC);

u = solve_lbvp(L,fn,B,BC.g,N); % Find Initial Condition for temperature

%% Plot Solution

U = reshape(u, [Grid.Ny, Grid.Nx, Grid.Nz]);

% Visualize the temperature distribution
figure;
slice(U, [], [], 1:Grid.Nz);
colorbar;
xlabel('Y-axis');
ylabel('X-axis');
zlabel('Z-axis');
title('Temperature Distribution in 3D');
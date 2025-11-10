function [Grid] = build_stokes_grid(Gridp)
% authors: Marc Hesse, Evan Carnahan

%% Build pressure grid

Gridp = build_grid_3D(Gridp)

%% Build x-velovity grid
Gridx.xmin = Gridp.xmin - 0.5*Gridp.dx; Gridx.xmax = Gridp.xmax + 0.5*Gridp.dx; Gridx.Nx = Gridp.Nx + 1;
Gridx.ymin = Gridp.ymin; Gridx.ymax = Gridp.ymax; Gridx.Ny = Gridp.Ny;
Gridx.zmin = Gridp.zmin; Gridx.zmax = Gridp.zmax; Gridx.Nz = Gridp.Nz;
Gridx = build_grid_3D(Gridx);

%% Build y-velocity grid
Gridy.xmin = Gridp.xmin; Gridy.xmax = Gridp.xmax; Gridy.Nx = Gridp.Nx;
Gridy.ymin = Gridp.ymin - 0.5*Gridp.dy; Gridy.ymax = Gridp.ymax + 0.5*Gridp.dy; Gridy.Ny = Gridp.Ny + 1;
Gridy.zmin = Gridp.zmin; Gridy.zmax = Gridp.zmax; Gridy.Nz = Gridp.Nz;
Gridy = build_grid_3D(Gridy);

%% Build z-velocity grid
Gridz.xmin = Gridp.xmin; Gridz.xmax = Gridp.xmax; Gridz.Nx = Gridp.Nx;
Gridz.ymin = Gridp.ymin; Gridz.ymax = Gridp.ymax; Gridz.Ny = Gridp.Ny;
Gridz.zmin = Gridp.zmin - 0.5*Gridp.dz; Gridz.zmax = Gridp.zmax + 0.5*Gridp.dz; Gridy.Nz = Gridp.Nz + 1;
Gridz = build_grid_3D(Gridz);

%% Define "Grid" Object that stores each grid

Grid.p = Gridp;
Grid.x = Gridx;
Grid.y = Gridy;
Grid.z = Gridz;

%% Helpful quantities
Grid.N = Grid.x.N + Grid.y.N + Grid.p.N + Grid.z.N; % total number of unknowns

%% Boundary dof's
% Unknown vector is ordered: u = [vx;vy;vz;p], but y-first with temperature?
% My build_grid_3D puts y-first so that's what I'm going to do, but let
% this be a record that before it was [vx;vy;p] for our 2D stokes problem

% Normal velocities on bnd's
Grid.dof_xmin_vx = Grid.x.dof_xmin + Grid.y.N;
Grid.dof_xmax_vx = Grid.x.dof_xmax + Grid.y.N;

Grid.dof_ymin_vy = Grid.y.dof_ymin;
Grid.dof_ymax_vy = Grid.y.dof_ymax;

Grid.dof_zmax_vz = Grid.z.dof_zmax + Grid.x.N + Grid.y.N;
Grid.dof_zmin_vz = Grid.z.dof_zmin + Grid.x.N + Grid.y.N;

% % Tangential velocities on bnd's (all)
Grid.dof_ymax_vx = Grid.x.dof_ymax + Grid.y.N;
Grid.dof_ymin_vx = Grid.x.dof_ymin + Grid.y.N;
Grid.dof_zmax_vx = Grid.x.dof_zmax + Grid.y.N;
Grid.dof_zmin_vx = Grid.x.dof_zmin + Grid.y.N;

Grid.dof_xmin_vy = Grid.y.dof_xmin;
Grid.dof_xmax_vy = Grid.y.dof_xmax;
Grid.dof_zmin_vy = Grid.y.dof_zmin;
Grid.dof_zmax_vy = Grid.y.dof_zmax;

Grid.dof_xmin_vz = Grid.z.dof_xmin + Grid.x.N + Grid.y.N;
Grid.dof_xmax_vz = Grid.z.dof_xmax + Grid.x.N + Grid.y.N;
Grid.dof_ymin_vz = Grid.z.dof_ymin + Grid.x.N + Grid.y.N;
Grid.dof_ymax_vz = Grid.z.dof_ymax + Grid.x.N + Grid.y.N;


% excluding extreme faces (THIS PROBABLY NEEDS TO BE FIXED)
Grid.dof_xmin_vt = [Grid.dof_xmin_vy(2:end-1); Grid.dof_xmin_vz(2:end-1)];
Grid.dof_xmax_vt = [Grid.dof_xmax_vy(2:end-1); Grid.dof_xmax_vz(2:end-1)];

Grid.dof_ymin_vt = [Grid.dof_ymin_vx(2:end-1); Grid.dof_ymin_vz(2:end-1)];
Grid.dof_ymax_vt = [Grid.dof_ymax_vx(2:end-1); Grid.dof_ymax_vz(2:end-1)];

Grid.dof_zmin_vt = [Grid.dof_zmin_vy(2:end-1); Grid.dof_zmin_vx(2:end-1)];
Grid.dof_zmax_vt = [Grid.dof_zmax_vy(2:end-1); Grid.dof_zmax_vx(2:end-1)];


% Pressures on bnd's
Grid.dof_xmin_p = Grid.p.dof_xmin + Grid.x.N + Grid.y.N + Grid.z.N;
Grid.dof_xmax_p = Grid.p.dof_xmax + Grid.x.N + Grid.y.N + Grid.z.N;
Grid.dof_ymin_p = Grid.p.dof_ymin + Grid.x.N + Grid.y.N + Grid.z.N;
Grid.dof_ymax_p = Grid.p.dof_ymax + Grid.x.N + Grid.y.N + Grid.z.N;
Grid.dof_zmin_p = Grid.p.dof_zmin + Grid.x.N + Grid.y.N + Grid.z.N;
Grid.dof_zmax_p = Grid.p.dof_zmax + Grid.x.N + Grid.y.N + Grid.z.N;

% Pressure constraint in center of domain
Grid.dof_pc = Grid.p.Nf + round(Grid.p.N/2);

%% Common useful BC's
% Penetration - set normal velocities on all boundaries to zero
Grid.dof_pene = [Grid.dof_xmin_vx; Grid.dof_xmax_vx; Grid.dof_ymin_vy; Grid.dof_ymax_vy; Grid.dof_zmin_vz; Grid.dof_zmax_vz];
Grid.N_pene = length(Grid.dof_pene);

% Slip - set all tangential velocities on all boundaries to zero
Grid.dof_slip = [Grid.dof_xmin_vt; Grid.dof_xmax_vt; Grid.dof_ymin_vt; Grid.dof_ymax_vt;Grid.dof_zmin_vt; Grid.dof_zmax_vt];
Grid.N_slip = length(Grid.dof_slip);

% Solid boundary - no slip and no penetration
Grid.dof_solid_bnd = unique([Grid.dof_pene;Grid.dof_slip]);
Grid.N_solid_bnd = length(Grid.dof_solid_bnd);

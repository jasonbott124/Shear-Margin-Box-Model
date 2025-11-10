function [D,Edot,Dp,Gp,Z,I,Mp,Gs] = build_stokes_ops(Grid) % repo
% authors: Marc Hesse, Evan Carnahan
% date: 4 Oct. 2019

%% Build operators for each variable
[Dp,Gp,~,Ip,Mp] = build_ops(Grid.p);
[Dx,Gx,~,Ix,~] = build_ops(Grid.x);
[Dy,Gy,~,Iy,~] = build_ops(Grid.y);

Nx = Grid.x.N;
Ny = Grid.y.N;
Nfx = Grid.x.Nfx;
Nfy = Grid.y.Nfy;

%% Extract x and y components of the velocity operators
Gxx = Gx(1:Nfx, :);
Gxy = Gx(Nfx+1:Grid.x.Nf, :);
Gyy = Gy(Grid.y.Nfx+1:Grid.y.Nf, :);
Gyx = Gy(1:Grid.y.Nfx, :);

Dxx = Dx(:, 1:Nfx);
Dxy = Dx(:, Nfx+1:Grid.x.Nf);
Dyy = Dy(:, Grid.y.Nfx+1:Grid.y.Nf);
Dyx = Dy(:, 1:Grid.y.Nfx);

% Zero blocks
Zxy = sparse(size(Gxx,1), size(Gyy,2));
Zyx = sparse(size(Gyy,1), size(Gxx,2));

%% Assemble Stokes operators
% Symmetric derivative
Edot = [Gxx,Zxy;...
        Zyx,Gyy;...
        0.5.*Gxy,0.5.*Gyx];

% Divergence of deviatoric stess tensor
D = [Dxx,Zyx',Dxy;...
Zxy',Dyy,Dyx];

Gs = -D';

% Zero block for the system matrix L
Z = sparse(size(Dp,1), size(Gp,2));

% Identity for all dof's
I = blkdiag(Ix, Iy, Ip);
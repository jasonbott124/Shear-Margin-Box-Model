function kappa_vec = set_kappa_field(u, Grid, Geom, kappa, kappa_sys)
    kappa_vec = nan(Grid.N,1);
    kappa_vec(Geom.seds.dof) = kappa.seds;
    kappa_vec(Geom.magma.dof) = kappa.magma;
    kappa_vec(Geom.basalt.dof) = kappa.basalt;
    kappa_vec(Geom.air.dof) = kappa.air;
    kappa_vec(Geom.ice.dof) = kappa_sys(u(Geom.ice.dof));
end
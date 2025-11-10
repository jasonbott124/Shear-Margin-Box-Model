function [phi,rho,cp,kappa,h,chi] = physical_properties_grl(Parameters)
% author: Marc A. Hesse
% date: 7 May 2018

% DESCRIPTION:
% Defines the physical properties of the phases in Ceres' volatile rich
% crust.

% INPUT:
% Ts = solidus temperature [K]
% Tl = liquidus temperature [K]
% L = latent heat of water [J/kg]
% phi = structure containing the initial volume fractions of phases [-]
%       phi.ice0, phi.hyd0, phi.sal0, phi.sil0

% OUTPUT:
% phi = structure containing the function handles for the volume fractions 
%       of the phases as function of temperature [-]
% rho = structure containing the constant densities of the phases [kg/m^3]
% cp = structure containing the function handles for the specific heat 
%      capacities of the phases as function of temperature [J/(kg K)]
% kappa = structure containing the function handles for the thermal 
%      conductivities of the phases as function of temperature [W/(m K)]
% h = structure containing the function handles for the specific enthalpy
%      of the phases as function of temperature [J/kg]
% chi = function handle of indicator function for melting interval [-]

% EXAMPLE CALL
% >> Ts = 245;
% >> Tl = 273;
% >> L = 3.34e5;
% >> phi0.ice = 0.3;
% >> phi0.hyd = 0.3;
% >> phi0.sal = 0.2;
% >> phi0.sil = 0.2;
% >> [phi,rho,cp,kappa,h,chi] = physical_properties_grl(Ts,Tl,L,phi0);

Tl = Parameters.Eutectic.T_liquidus; 
Ts = Parameters.Eutectic.T_solidus; 
Lw = Parameters.Lw;
phi = Parameters.phi;

%% Indicator function for the melting interval [-]
% This is simply a function that is 1 for Ts < T < Tl and 0 elsewhere
chi.ice_water = @(T) (T > Ts).*(T < Tl);

% chi.magma_basalt = @(T) (T > Ts.basalt_magma).*(T < Tl.basalt_magma);
% chi.seds = @(T) (T>Ts.seds).*(T<Tl.seds);
% chi.air = @(T) (T>Ts.air).*(T<Tl.air);

%% Densities of the phases [kg/m^3]
rho.ice   = 917;   % Density of Ice [kg/m^3]
rho.wat = 1e3;   % Density of Water [kg/m^3]

% rho.magma = 2.4e3; % Density of Si-Rich Trachyte Magma [kg/m^3]
% rho.seds = 2.2e3; % Density of "Sediments"
% rho.basalt = 2.9e3; % Density of Basalt (According to Google)
% rho.air = 1.225;

%% Volume fractions of the phases as function of T [-]

f = @(T) (Tl-T)/(Tl-Ts); % function decreasing linearly from 1 to 0 between Ts and Tl

phi.wat_star = phi.ice0; % what is this???
phi.ice = @(T) phi.ice0*(T<=Ts) + phi.ice0*f(T).*chi.ice_water(T);
phi.wat = @(T) phi.wat_star*(T >= Tl) + phi.wat_star*(1-f(T)).*chi.ice_water(T);

% phi.basalt = @(T) phi.basalt0*(T<=Ts) + phi.basalt0*f(T).*chi.basalt_magma(T);
% Setting melt fractions of the rock equal to zero (should it be zero or 1?)
% phi.seds = 0;
% phi.air = 0;
% phi.magma = @(T) phi.magma_star*(T >= Tl.basalt_magma)+phi.magma_star*(1-f(T)).*chi.basalt_magma(T);

%% Heat capacity [J/kg/K]
% All phases considered here have at most a linear variation in cp
heat_capacity = @(T,const) const.a + const.b*T;

% Coefficients in the cp polynomial
cp.const.ice.a = 185;  cp.const.ice.b = 7.037;
cp.const.wat.a = 4200; cp.const.wat.b = 0;

% cp.const.sil.a = 2000; cp.const.sil.b = 0;
% cp.const.sal.a = 920;  cp.const.sal.b = 0;
% cp.const.hyd.a = 494;  cp.const.hyd.b = 6.1;

% cp.const.seds.a = 2000; cp.const.seds.b = 0; 
% cp.const.magma.a = 850; cp.const.magma.b = 0;
% cp.const.basalt.a = 1000; cp.const.basalt.b = 0;
% cp.const.air.a = 1005; cp.const.air.b = 0;

% Variable Heat Capacity
cp.ice = @(T) heat_capacity(T,cp.const.ice);
cp.wat = @(T) heat_capacity(T,cp.const.wat);
% cp.seds = @(T) heat_capacity(T,cp.const.seds);
% cp.basalt = @(T) heat_capacity(T,cp.const.basalt);
% cp.magma = @(T) heat_capacity(T,cp.const.magma);
% cp.air = @(T) heat_capacity(T,cp.const.air);

%% Thermal Conductivity [W/m/K]
% kappa.basalt = 2.65 ;  % Thermal Conductivity of Basalt [W/m K]
kappa.ice = Parameters.kappa.ice_wolf; 
kappa.wat = @(T) 0.56 + T*0;
% kappa.seds =   2.8 ;  % Thermal Conductivity of Sediment [W/m K]
% kappa.magma =   1.6 ;   % Thermal Conductivity of Magma [W/m K]
% kappa.air =   3.0 ;   % Thermal Conductivity of Air [W/m K]

%% Specific enthalpy of the phases
% General integral of the linear heat capacity function
integral_heat_capacity = @(T,const,material,Ts) const.a*(T - Ts) + const.b/2 * (T.^2 - Ts^2);

h.ice = @(T) integral_heat_capacity(T,cp.const.ice,Ts);
h.wat = @(T) integral_heat_capacity(T,cp.const.wat) + Lw;

% h.air = @(T) integral_heat_capacity(T,cp.const.air,Ts.air);
% h.basalt = @(T) integral_heat_capacity(T,cp.const.basalt,Ts.basalt_magma);
% h.seds = @(T) integral_heat_capacity(T,cp.const.seds,Ts.seds);
% h.magma = @(T) integral_heat_capacity(T,cp.const.magma,Ts.basalt_magma) + Lm;


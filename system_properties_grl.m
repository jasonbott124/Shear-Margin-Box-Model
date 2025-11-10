function [phi,kappa_sys,rhocp_sys,H,rho,cp,kappa,h,chi] = system_properties_grl(Parameters)
% author: Marc A. Hesse
% date: 8 May 2018

% DESCRIPTION:
% Defines the physical properties of five phase system composing Ceres'
% crust.

% INPUT:
% Ts = solidus temperature [K]
% Tl = liquidus temperature [K]
% L = latent heat of water [J/kg]

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
    
Tl = Parameters.Eutectic.T_liquidus; 
Ts = Parameters.Eutectic.T_solidus;

[phi,rho,cp,kappa,h,chi] = physical_properties_grl(Parameters);

DT = Tl - Ts;  % Melting interval

%% Total enthalpy of the system

H = @(T) rho.ice*phi.ice(T).*h.ice(T) + ...
         rho.wat*phi.wat(T).*h.wat(T) ;

% H_bm_sys = @(T) rho.basalt*phi.basalt(T).*h.basalt(T) + ...
%                 rho.magma*phi.magma(T).*h.magma(T);

         % rho.magma*phi.magma.*h.magma(T) + ...
         % rho.seds*phi.seds.*h.seds(T) + ...
         % rho.basalt*phi.basalt.*h.basalt(T) + ...
         % rho.air*phi.air.*h.air(T) ;

kappa_sys = @(T) phi.ice(T).*kappa.ice(T) + ...
                 phi.wat(T).*kappa.wat(T);

% kappa_bm_sys = @(T) phi.basalt(T).*kappa.basalt(T) + ...
%                     phi.magma(T).*kappa.magma(T);
                 
                 % phi.seds.*kappa.seds(T) + ...
                 % phi.magma.*kappa.magma(T) + ...
                 % phi.basalt.*kappa.basalt(T) + ...
                 % phi.air.*kappa.air(T) + ...

%% Effective heat capacity of system
rhocp_mean = @(T) rho.ice*phi.ice(T).*cp.ice(T) + ...
                  rho.wat*phi.wat(T).*cp.wat(T);

% rhocp_bm_mean = @(T) rho.basalt*phi.basalt(T).*cp.basalt(T) + ...
%                   rho.magma*phi.magma(T).*cp.magma(T);

                  % rho.seds*phi.seds.*cp.seds(T) + ...
                  % rho.magma*phi.magma.*cp.magma(T) + ...
                  % rho.basalt*phi.basalt.*cp.basalt(T) + ...
                  % rho.basalt*phi.air.*cp.air(T) ;


latent_heat = @(T) chi(T)/DT.*(rho.wat*phi.wat_star*h.wat(T)...
                              -rho.ice*phi.ice0*h.ice(T));

% latent_heat = @(T) chi(T)/DT.*(rho.wat*phi.wat_star*h.wat(T)...
%                               -rho.ice*phi.ice0*h.ice(T));
                              
                              % -rho.seds*phi.seds0*h.seds(T)...
                              % -rho.magma*phi.magma0*h.magma(T)...
                              % -rho.basalt*phi.basalt0*h.basalt(T) ...
                              % -rho.air*phi.air0*h.air(T)); 
 
                              % should latent heat include the other materials?                             


rhocp_sys = @(T) rhocp_mean(T) + latent_heat(T);


end
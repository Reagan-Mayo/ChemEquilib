%% main.m
% 11/17/2022
clear; clc; close all;
%% PURPOSE
% Main code to collect JANAF data, iterate until converge at a chemical
% equilibrium, and output figures.
%% INPUTS

% Species Information
spec = {'H2','O2','N2','H2O','OH','O','H','NO','Ne'}; %species under consideration
spec_Kp = {'H2O','OH','NO','H','O'};   %species that require their Kp
X.Ne = 0.01;    %Neon mole fraction = 0

% Inlet Conditions
f.p = 20;   %fuel pressure [atm]
f.T = 298;  %fuel temp [K]
o.p = f.p;   %oxidizer pressure [atm]
o.T = 500;  %oxidizer temp [K]

% Combustion Properties
phi = 0.7:0.1:1.3;

% Convergence Criteria
eps = 1e-6;

%% EXECUTE

for j = 1:length(phi)
    % Collect JANAF thermochemical data
    for i = 1:length(spec_Kp)
        [T.(spec_Kp{i}), Kp.(spec_Kp{i})] = readJANAF(spec_Kp{i});
    end
    
    % Iterate equilibrium problem until converged
    
    err = inf;
    Tguess = 2000;  %guessed adiabatic flame temperature [K]
    while err >= eps
        [Te,pe,X] = thermoChemEquilib(Tguess,f.p,T,Kp,spec_Kp,phi(j),X.Ne);
        err = abs(Te - Tguess);
        Tguess = Te;
    end
    
    % Save data
    T_save(j) = Te; %adiabatic flame temp [K]
    X_save(j) = X;  %Mole fractions
end

% Plot




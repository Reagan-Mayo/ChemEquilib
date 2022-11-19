%% main.m
% 11/17/2022
clear; clc; close all;
%% PURPOSE
% Main code to collect JANAF data, iterate until converge at a chemical
% equilibrium, and output figures.
%% INPUTS

% Species Information
spec = {'H2','O2','N2','H2O','OH','O','H','NO','Ne'}; %species under consideration
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

% Convergence Rate
cr = 0.2;

% Adiabatic Flame Temperature Initial Guess
Tguess = 2000;  %[K]

%% EXECUTE

for j = 1:length(phi)
    % Collect JANAF thermochemical data
    for i = 1:length(spec)
        [T.(spec{i}), h_hTref.(spec{i}), dhf.(spec{i}), Kp.(spec{i})] = readJANAF(spec{i});
    end
    
    % Iterate equilibrium problem until converged
    % Prepare Loop
    Tguess1 = inf;
    err = inf; err1 = inf;
    cr1 = cr;
    while err >= eps    % While convergence not met
        [X,err] = thermoChemEquilib(f,o,Tguess,f.p,T,h_hTref,dhf,Kp,spec,phi(j),X.Ne);
        
        if err*err1 < 0     %if error passes 0, change direction of temperature steps
            cr1 = -0.1*cr1;
        end
        
        % Prepare for next loop
        err1 = err;
        Tguess1 = Tguess;
        Tguess = Tguess*(1+cr1);
    end
    
    % Outlet Conditions
    pe = o.p;
    Te = Tguess;

    % Save data
    T_save(j) = Te; %adiabatic flame temp [K]
    X_save(j) = X;  %Mole fractions
end

% Plot




%% readJANAF.m
% 11/17/2022
%% Purpose
% Read local JANAF data tables for a given species.
%% I/O
% INPUT
%   - spec: string of species
% OUTPUT
%   - T: temperature [K]
%   - Cp: specific heat at constant pressure [K/(K mol)]
%   - Kf: equilibrium constant of formation [-]

%% EXECUTE

function [T, Kf] = readJANAF(spec)

fileStr = [spec,'.txt'];

data = readmatrix(fileStr);

T = data(:,1);
logKf = data(:,end);
% Kf = exp(logKf);   %TODO: ensure log is ln() and not log10()
Kf = 10.^(logKf);
end
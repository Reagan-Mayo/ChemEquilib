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

function [T, h_hTref, dhf, Kf] = readJANAF(spec)

fileStr = [spec,'.txt'];

data = readmatrix(fileStr);

T = data(:,1);          %[K]
h_hTref = data(:,5);    %[kJ/mol]
dhf = data(:,6);        %[kJ/mol]
logKf = data(:,end);
Kf = 10.^(logKf);    %TODO: see if log is ln or log10
end
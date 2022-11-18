%% thermoChemEquilib.m
% 11/17/2022
%% PURPOSE
% to perform single pass calculation of thermochemical equilibrium.
%% I/O
% INPUT
%   - T: (guessed) equilibrium temperature [K]
%   - p: equilibrium pressure [atm]
%   - Tvec: temperature vector from JANAF table for all relevant species
%   - Kpvec: Kp vector from JANAF table for all relevant species
%   - spec_Kp: species that requrie Kp values
%   - phi: equivalence ratio
%   - XNe: Neon mole fraction ratio
% OUPUT
%   - Te: new equilibrium temperature [K]
%   - pe: equilibrium pressure [atm]
%   - X: equilibrium mole fractions
%% EXECUTE

function [Te,pe,X] = thermoChemEquilib(T,p,Tvec,Kpvec,spec_Kp,phi,XNe)

% Interpolate Kp's
for i = 1:length(spec_Kp)
    Kp.(spec_Kp{i}) = interp1(Tvec.(spec_Kp{i}),Kpvec.(spec_Kp{i}),T);
end
Kp.Ne = 1;  %inert

% Formation Reaction Equations

syms XH2O XH2 XO2 XOH XNO XN2 XH XO

eqn1 = Kp.H2O == XH2O./(XH2.*(XO2.*p).^(1/2));
eqn2 = Kp.OH == XOH ./ (XH2.*XO2).^(1/2);
eqn3 = Kp.NO == XNO ./ (XN2.*XO2).^(1/2);
eqn4 = Kp.H == XH ./ (XH2.*p).^(1/2);
eqn5 = Kp.Ne == 1;

% Atom Conservation Equations
eqn6 = 2.*phi == (2.*XH2 + 2.*XH2O + XOH + XH)./(2.*XO2 + XH2O + XOH + XO + XNO);
eqn7 = 3.76 == (2.*XN2 + XNO)./(2.*XO2 + XH2O + XOH + XO + XNO);
eqn9 = 1 == XH2 + XO2 + XN2 + XH2O + XOH + XO + XH + XNO + XNe;
if XNe ~= 0
    eqn8 = XNe == XNe ./ (2.*XO2 + XH2O + XOH + XO + XNO);
else
    eqn8 = Kp.O == XO ./ (XO2.*p).^(1/2);
end

S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn6,eqn7,eqn8,eqn9],[XH2O,XH2,XO2,XOH,XNO,XN2,XH,XO],[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]) %TODO: FIX THIS, NOT SOLVING

end
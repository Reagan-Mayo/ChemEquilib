%% thermoChemEquilib.m
% 11/17/2022
%% PURPOSE
% to perform single pass calculation of thermochemical equilibrium.
%% I/O
% INPUT
%   - f: fuel inlet conditions
%   - o: oxidizer inlet conditions
%   - T2: (guessed) equilibrium temperature [K]
%   - p2: equilibrium pressure [atm]
%   - Tvec: temperature vector from JANAF table for all relevant species
%   - Kpvec: Kp vector from JANAF table for all relevant species
%   - spec: species that requrie Kp values
%   - phi: equivalence ratio
%   - XNe: Neon mole fraction ratio
% OUPUT
%   - Te: new equilibrium temperature [K]
%   - pe: equilibrium pressure [atm]
%   - X: equilibrium mole fractions
%   - TadEqnIneq: Inequality (rxnt - prod) in the total enthalpy equation
%% EXECUTE

function [X,TadEqnIneq] = thermoChemEquilib(f,o,T2,p2,Tvec,h_hTrefvec,dhfvec,Kpvec,spec,phi,XNe)

% Interpolate Kp's
for i = 1:length(spec)
    Kp.(spec{i}) = interp1(Tvec.(spec{i}),Kpvec.(spec{i}),T2);
end
Kp.Ne = 1;  %inert

% Formation Reaction Equations

syms XH2O XH2 XO2 XOH XNO XN2 XH XO

eqn1 = log(Kp.H2O) == log(XH2O./(XH2.*(XO2.*p2).^(1/2)));
eqn2 = log(Kp.OH) == log(XOH ./ (XH2.*XO2).^(1/2));
eqn3 = log(Kp.NO) == log(XNO ./ (XN2.*XO2).^(1/2));
eqn4 = log(Kp.H) == log(XH ./ (XH2.*p2).^(1/2));
eqn5 = log(Kp.Ne) == log(1);

% Atom Conservation Equations
eqn6 = 2.*phi == (2.*XH2 + 2.*XH2O + XOH + XH)./(2.*XO2 + XH2O + XOH + XO + XNO);
eqn7 = 3.76 == (2.*XN2 + XNO)./(2.*XO2 + XH2O + XOH + XO + XNO);
eqn9 = 1 == XH2 + XO2 + XN2 + XH2O + XOH + XO + XH + XNO + XNe;

if XNe ~= 0
    eqn8 = XNe == XNe ./ (2.*XO2 + XH2O + XOH + XO + XNO);
else
    eqn8 = Kp.O == XO ./ (XO2.*p2).^(1/2);
%     eqn8 = Kp.O == abs(XO) ./ (abs(XO2).*p2).^(1/2);
end

S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn6,eqn7,eqn8,eqn9],[XH2O,XH2,XO2,XOH,XNO,XN2,XH,XO],[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]);

% Store Mole Fraction Ratio Solutions
fn = fieldnames(S);
for i = 1:length(fn)
    X.(fn{i}(2:end)) = double(S.(fn{i}));
end
X.Ne = XNe;

% Interpolate Products Enthalpies at Guessed Equilibrium Temperatures
for i = 1:length(spec)
    h_hTref.(spec{i}) = interp1(Tvec.(spec{i}),h_hTrefvec.(spec{i}),T2);
    dhf.(spec{i}) = interp1(Tvec.(spec{i}),dhfvec.(spec{i}),T2);
    prod.(spec{i}) = X.(spec{i}) .* (h_hTref.(spec{i}) + dhf.(spec{i}));
end

% Interpolate Reactants Enthalpies at Inlet Temperatures
spec_ox = {'O2','N2','Ne'};
spec_f = {'H2'};

for i = 1:length(spec_ox)   %Oxidizer and related gases
    h_hTrefRxnt.(spec_ox{i}) = interp1(Tvec.(spec_ox{i}),h_hTrefvec.(spec_ox{i}),o.T);
    dhfRxnt.(spec_ox{i}) = interp1(Tvec.(spec_ox{i}),dhfvec.(spec_ox{i}),o.T);
    rxnt.(spec_ox{i}) = (h_hTrefRxnt.(spec_ox{i}) + dhfRxnt.(spec_ox{i}));
end
for i = 1:length(spec_f)    %Fuel
    h_hTrefRxnt.(spec_f{i}) = interp1(Tvec.(spec_f{i}),h_hTrefvec.(spec_f{i}),f.T);
    dhfRxnt.(spec_f{i}) = interp1(Tvec.(spec_f{i}),dhfvec.(spec_f{i}),f.T);
    rxnt.(spec_f{i}) = (h_hTrefRxnt.(spec_f{i}) + dhfRxnt.(spec_f{i}));
end

% Solve for Inequality
LHS = phi.*rxnt.H2 + 0.5.*(rxnt.O2 + 3.76.*rxnt.N2) + X.Ne .* rxnt.N2;
RHS = prod.H2 + prod.O2 + prod.N2 + prod.H2O + prod.OH + prod.O + prod.H + prod.NO + prod.Ne;
TadEqnIneq = LHS - RHS;



end
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
format long
%% New Way
% X = ['H2','O2','N2','H2O','OH','O','H','NO'];
% X0 = [1 1 1 1 1 1 1 1];
% for i = 1:length(spec)-1
%     c(i) = Kp.(spec{i});
% end
% c(9:11) = [phi p2 XNe];
% fun = @(X)chemCompSolver(X,c)
% 
% options = optimoptions(@fsolve,'MaxFunctionEvaluations',2e4)
% S = fsolve(fun,X0,options)



%% Old Way
syms XH2O XH2 XO2 XOH XNO XN2 XH XO

eqn1 = log10(Kp.H2O) == log10(XH2O./(XH2.*(XO2.*p2).^(1/2)));
eqn2 = log10(Kp.OH) == log10(XOH ./ (XH2.*XO2).^(1/2));
eqn3 = log10(Kp.NO) == log10(XNO ./ (XN2.*XO2).^(1/2));
eqn4 = log10(Kp.H) == log10(XH.*p2.^0.5 ./ (XH2).^(1/2));
eqn5 = log10(Kp.Ne) == log10(1);

% Atom Conservation Equations
eqn6 = 2.*phi == (2.*XH2 + 2.*XH2O + XOH + XH)./(2.*XO2 + XH2O + XOH + XO + XNO);
eqn7 = 3.76 == (2.*XN2 + XNO)./(2.*XO2 + XH2O + XOH + XO + XNO);
eqn9 = 1 == XH2 + XO2 + XN2 + XH2O + XOH + XO + XH + XNO + XNe;
eqn8 = log10(Kp.O) == log10(XO.*p2.^0.5 ./ (XO2).^(1/2));



assume(XH2O,{'real','positive'}); assume(XH2,{'real','positive'}); assume(XO2,{'real','positive'}); assume(XOH,{'real','positive'});
assume(XNO,{'real','positive'}); assume(XN2,{'real','positive'}); assume(XH,{'real','positive'}); assume(XO,{'real','positive'});

S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn6,eqn7,eqn8,eqn9],[XH2O,XH2,XO2,XOH,XNO,XN2,XH,XO],[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]);
if isempty(S.XH2)
    S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn6,eqn7,eqn8,eqn9],[XH2O,XH2,XO2,XOH,XNO,XN2,XH,XO]);
end
%% End Old Way


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
ntot = 3.76/(2*X.N2 + X.NO);
LHS = phi.*rxnt.H2 + 0.5.*(rxnt.O2 + 3.76.*rxnt.N2) + X.Ne .* rxnt.N2;
RHS = ntot*prod.H2 + ntot*prod.O2 + ntot*prod.N2 + ntot*prod.H2O + ntot*prod.OH + ntot*prod.O + ntot*prod.H + ntot*prod.NO + ntot*prod.Ne;
TadEqnIneq = LHS - RHS;



end
%% plotThermEquilib.m
% 11/19/2022
%% PURPOSE
% 
%% INPUT

Ts1 = load("Ts1.mat");
Ts2 = load("Ts2.mat");
Ts3 = load("Ts3.mat");
Ts4 = load("Ts4.mat");

Xs1 = load("Xs1.mat");
Xs2 = load("Xs2.mat");
Xs3 = load("Xs3.mat");
Xs4 = load("Xs4.mat");

spec = {'H2','O2','N2','H2O','OH','O','H','NO','Ne'}; %species under consideration
phi = 0.7:0.1:1.3;

%% EXECUTE

% Scenario 1
figure
for i = 1:length(spec)
    chi = [];
    for j = 1:length(phi)
        chi(j) = Xs1.X_save(j).(spec{i});
    end
    semilogy(phi,chi,'LineWidth',1); hold on;
end
legend(spec)
xlabel('Equivalence Ratio, \phi')
ylabel('Mole Fraction, \chi')
xlim([0.7,1.3])
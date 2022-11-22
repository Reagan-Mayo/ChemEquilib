%% plotThermEquilib.m
% 11/19/2022
%% PURPOSE
% 
%% INPUT

clear; clc; close all;

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

%% X vs phi @ p = 20atm
figure
for i = 1:length(spec)
    chi = [];
    for j = 1:length(phi)
        chi(j) = Xs1.X_save(j).(spec{i});
    end
    semilogy(phi,chi,'LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
end
legend(spec,'Location', 'bestoutside')
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Equivalence Ratio, \phi','FontSize',14)
ylabel('Mole Fraction, \chi','FontSize',14)
xlim([0.7,1.3])
grid on

hold off

%% T vs phi @ p = 20atm
figure
plot(phi,Ts1.T_save,'LineWidth',2,'Marker','d','MarkerFaceColor','auto')
xlabel('Equivalence Ratio, \phi')
ylabel('Adiabatic Flame Temperature [K]')
xlim([0.7,1.3])
grid on
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
hold off

%% T vs p @ phi = 1 @ phi = 1
figure
p = [2 10 20]; T = [Ts3.T_save, Ts2.T_save, Ts1.T_save(4)];
plot(p,T,'LineWidth',2,'Marker','d','MarkerFaceColor','auto')
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Pressure [atm]')
ylabel('Adiabatic Flame Temperature [K]')
xlim([0 22])
grid on
hold off

%% NOx vs phi @ p = 20atm
figure
for i = 1:length(phi)
    chi(i) = Xs1.X_save(i).NO;
end
plot(phi,chi,'LineWidth',2,'Marker','d','MarkerFaceColor','auto')
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Equivalence Ratio, \phi')
ylabel('NOx Concentration, \chi_{NO}')
xlim([0.7 1.3])
grid on
hold off

%% NOx vs p @ phi = 1
figure
p = [2 10 20];
chi = [Xs3.X_save.NO, Xs2.X_save.NO, Xs1.X_save(4).NO];
plot(p,chi,'LineWidth',2,'Marker','d','MarkerFaceColor','auto')
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Pressure [atm]')
ylabel('NOx Concentration, \chi_{NO}')
xlim([0 22])
ylim([3e-3, 4e-3])
grid on
hold off

%% phi vs X, Code vs GasEQ
geq{1} = [0.68559 0.25489 0.05172 1.95e-4 2.11e-3 1.24e-5 9.27e-5 0.00539];
geq{2} = [0.67238 0.2846 0.03225 7.42e-4 3.75e-3 5.93e-5 2e-4 0.00602];
geq{3} = [0.65979 0.31159 0.01482 2.61e-3 0.00533 2.25e-4 2.99e-4 0.00534];
geq{4} = [0.64676 0.33097 2.79e-3 0.01072 0.00505 7.28e-4 2.19e-4 2.76e-3];
geq{5} = [0.62879 0.33029 2.59e-4 0.0356 2.81e-3 1.34e-3 6.78e-5 8.33e-4];
geq{6} = [0.60922 0.32188 5.04e-5 0.06531 1.64e-3 1.54e-3 2.48e-5 3.39e-4];
geq{7} = [0.59041 0.31269 1.491e-5 0.09415 1.03e-3 1.52e-3 1.08e-5 1.68e-4];
spec_geq = {'N2','H2O','O2','H2','OH','H','O','NO'};
for i = 1:7
    for j = 1:length(spec_geq)
        Xgeq(i).(spec_geq{j}) = geq{i}(j);
    end
end

figure
color = {"#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#A2142F", "#0000FF", "#000000"};
for j = 1:length(spec_geq)
    for i = 1:7
        chi(i) = Xs1.X_save(i).(spec_geq{j});
        chigeq(i) = Xgeq(i).(spec_geq{j});
    end
    p(j) = semilogy(phi,chi,'-','Color',color{j},'LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
    q(j) = semilogy(phi,chigeq,'--','Color',color{j},'LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
end
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Equivalence Ratio, \phi')
ylabel('Concentration, \chi')
legend([p(1:8)],[spec_geq],'Location', 'bestoutside')
xlim([0.7 1.3])
grid on
hold off

%% phi vs T, Code vs GasEQ

Tgeq = [2141.1 2302.3 2446.6 2553.2 2556.3 2517.0 2472.6];
figure
plot(phi,Ts1.T_save,'-','LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
plot(phi,Tgeq,'--','LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Equivalence Ratio, \phi')
ylabel('Adiabatic Flame Temperature [K]')
legend({'Present Code','GasEQ'})
xlim([0.7 1.3])
grid on
hold off

%% p vs X, Code vs GasEQ @ phi = 1

pressure = [2 10 20];
x = [Xs3.X_save Xs2.X_save Xs1.X_save(4)];
geq{1} = [0.64292 0.32009 0.00514 0.01735 0.00833 2.17e-3 6.73e-4 3.32e-3];
geq{2} = [0.64575 0.32808 3.4e-3 0.01252 0.00595 1.03e-3 3.14e-4 2.95e-3];
geq{3} = [0.64676 0.33097 2.79e-3 0.01072 0.00505 7.28e-4 2.19e-4 2.76e-3];
Xgeq = [];
for j = 1:length(spec_geq)
    X.(spec_geq{j}) = [x(1).(spec_geq{j}) x(2).(spec_geq{j}) x(3).(spec_geq{j})];
    Xgeq.(spec_geq{j}) = [geq{1}(j) geq{2}(j) geq{3}(j)];
end

figure
for j = 1:length(spec_geq)
    p(j) = semilogy(pressure,X.(spec_geq{j}),'-','Color',color{j},'LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
    q(j) = semilogy(pressure,Xgeq.(spec_geq{j}),'--','Color',color{j},'LineWidth',2,'Marker','d','MarkerFaceColor','auto'); hold on;
end
ax = gca;
ax.FontSize = 13;
ax.FontName = 'Times New Roman';
xlabel('Pressure [atm]')
ylabel('Concentration, \chi')
legend([p(1:8)],[spec_geq],'Location', 'bestoutside')
xlim([0 22])
grid on
hold off



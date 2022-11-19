function [F] = chemCompSolver(x,c)

spec = {'H2','O2','N2','H2O','OH','O','H','NO'}; %species under consideration
for i = 1:length(spec)
    Kp.(spec{i}) = c(i);
end

XH2 = x(1);
XO2 = x(2);
XN2 = x(3);
XH2O = x(4);
XOH = x(5);
XO = x(6);
XH = x(7);
XNO = x(8);

phi = c(9);
p2 = c(10);
XNe = c(11);

if XNe == 0
    F = [log10(Kp.H2O) - log10(XH2O./(XH2.*(XO2.*p2).^(1/2)));
        log10(Kp.OH) - log10(XOH ./ (XH2.*XO2).^(1/2));
        log10(Kp.NO) - log10(XNO ./ (XN2.*XO2).^(1/2));
        log10(Kp.H) - log10(XH .* p2.^(0.5) ./ (XH2).^(1/2));
        2.*phi - (2.*XH2 + 2.*XH2O + XOH + XH)./(2.*XO2 + XH2O + XOH + XO + XNO);
        3.76 - (2.*XN2 + XNO)./(2.*XO2 + XH2O + XOH + XO + XNO);
        1 - XH2 + XO2 + XN2 + XH2O + XOH + XO + XH + XNO + XNe;
        log10(Kp.O) - log10(XO.*p2.^0.5 ./ (XO2).^(1/2))];
else
    F = [log10(Kp.H2O) - log10(XH2O./(XH2.*(XO2.*p2).^(1/2)));
        log10(Kp.OH) - log10(XOH ./ (XH2.*XO2).^(1/2));
        log10(Kp.NO) - log10(XNO ./ (XN2.*XO2).^(1/2));
        log10(Kp.H) - log10(XH .*p2^0.5./ (XH2).^(1/2));
        2.*phi - (2.*XH2 + 2.*XH2O + XOH + XH)./(2.*XO2 + XH2O + XOH + XO + XNO);
        3.76 - (2.*XN2 + XNO)./(2.*XO2 + XH2O + XOH + XO + XNO);
        1 - XH2 + XO2 + XN2 + XH2O + XOH + XO + XH + XNO + XNe;
        XNe - XNe ./ (2.*XO2 + XH2O + XOH + XO + XNO)];
end


end


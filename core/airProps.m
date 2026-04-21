

function [Kaire, den_Aire, Cp_aire, Vd_Aire, Vs_Aire, Alfa_aire, Pr_Aire, BET_A] = airProps(Taire_K)
% airProps  Propiedades del aire seco en función de la temperatura
%
% INPUT
%   Taire_K : temperatura absoluta del aire [K]
%
% OUTPUT
%   Kaire    : conductividad térmica [W/m·K]
%   den_Aire : densidad [kg/m³]
%   Cp_aire  : calor específico a presión constante [J/kg·K]
%   Vd_Aire  : viscosidad dinámica [Pa·s]
%   Vs_Aire  : viscosidad cinemática [m²/s]
%   Alfa_aire: difusividad térmica alpha = k/(rho*Cp) [m²/s]
%   Pr_Aire  : número de Prandtl (-)
%   BET_A    : coeficiente de expansión térmica beta = 1/T [1/K]

    % pasar de K a °C para las correlaciones
    Tpe = Taire_K ;        % [°C]

    % conductividad térmica
    Kaire = 0.0244 + (0.6763e-4).*Tpe;

    % densidad (P/(R*T)), reescrito como 353.44/(T[°C]+273.15)
    den_Aire = 353.44./(Tpe + 273.15);   % [kg/m³]

    % Cp aire (correlación polinómica)
    Cp_aire = 999.2 + 0.1434.*Tpe + (1.101e-4).*Tpe.^2;   % [J/kg·K]

    % viscosidad dinámica
    Vd_Aire = (1.719e-5) + (4.620e-8).*Tpe;               % [Pa·s]

    % viscosidad cinemática
    Vs_Aire = Vd_Aire ./ den_Aire;                        % [m²/s]

    % difusividad térmica
    Alfa_aire = Kaire ./ (den_Aire .* Cp_aire);           % [m²/s]

    % número de Prandtl
    Pr_Aire = Vs_Aire ./ Alfa_aire;                       % [-]

    % coeficiente de expansión térmica
    BET_A = 1 ./ Taire_K;                                 % [1/K]
end
function props = humidAirProps_T_RH(Tdb_C, RH_frac, Patm_kPa)
% humidAirProps_T_RH  Propiedades de aire húmedo sin tablas
%
% INPUTS
%   Tdb_C     : temperatura de bulbo seco [°C]
%   RH_frac   : humedad relativa en fracción (0–1)
%   Patm_kPa  : presión atmosférica [kPa]
%
% OUTPUT (struct)
%   props.w    : humedad específica [kg_v/kg_aire seco]
%   props.h    : entalpía específica [kJ/kg_aire seco]
%   props.v    : volumen específico [m³/kg_aire seco]
%   props.Pw   : presión parcial de vapor [kPa]
%   props.Pa   : presión parcial de aire seco [kPa]
%   props.Tdp  : temperatura de rocío [°C]
%   props.Pws  : presión de saturación a Tdb [kPa]

    % Constantes
    Ra = 0.287042;  % m³·kPa/(kg·K)
    
    % 1) Presión de saturación del agua a Tdb (Magnus/Tetens)
    %    válido típicamente 0–60 °C
    Pws = 0.61078 .* exp(17.27 .* Tdb_C ./ (Tdb_C + 237.3));  % [kPa]
    
    % 2) Presión de vapor y de aire seco
    Pw = RH_frac .* Pws;          % [kPa]
    Pa = Patm_kPa - Pw;           % [kPa]
    
    % 3) Humedad específica (kg_vapor / kg_aire seco)
    w  = 0.622 .* (Pw ./ Pa);
    
    % 4) Entalpía (kJ/kg aire seco) (ASHRAE aprox)
    %    Cp_air ~ 1.005 kJ/kgK ; h_fg(T) ~ 2501 + 1.86*T
    h = 1.005 .* Tdb_C + w .* (2501 + 1.86 .* Tdb_C);
    
    % 5) Volumen específico (m³/kg aire seco)
    T_K = Tdb_C + 273.15;
    v   = Ra .* T_K ./ Pa;
    
    % 6) Punto de rocío (inversión de Magnus/Tetens usando Pw)
    %    Pw = 0.61078 * exp(17.27*Tdp/(Tdp+237.3))
    gamma = log(Pw ./ 0.61078);
    Tdp   = 237.3 .* gamma ./ (17.27 - gamma);
    
    % Empaquetar
    props.w   = w;
    props.h   = h;
    props.v   = v;
    props.Pw  = Pw;
    props.Pa  = Pa;
    props.Tdp = Tdp;
    props.Pws = Pws;
end
function [cp_mass, cp_molar] = Cp_h2o(T_C)
    % --- constants ---
    R = 8.31446261815324;     
    M = 18.01528e-3;          

    T = T_C(:) + 273.15;      % K

      a_low  = [ ...
        4.19864056; ...
       -2.03643410e-3; ...
        6.52040211e-6; ...
       -5.48797062e-9; ...
        1.77197817e-12; ...
       -3.02937267e4; ...
       -0.849032208 ];

    a_high = [ ...
        3.03399249; ...
        2.17691804e-3; ...
       -1.64072518e-7; ...
       -9.70419870e-11; ...
        1.68200992e-14; ...
       -3.00042971e4; ...
        4.96677010 ];

    % --- sanity checks ---
    if any(~isfinite(T))
        error('Non-finite temperature detected.');
    end
    if any(T < 200 | T > 3500)
        warning('Some T values are outside NASA7 validity range (200–3500 K). Results may be invalid.');
    end

    % --- piecewise evaluation Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4 ---
    cpR = zeros(size(T));

    idxLow  = (T < 1000);     % includes your 90–400°C range (363–673 K)
    idxHigh = ~idxLow;

    cpR(idxLow)  = polyCpR(T(idxLow),  a_low);
    cpR(idxHigh) = polyCpR(T(idxHigh), a_high);

    cp_molar = cpR * R;       % J/(mol·K)
    cp_mass  = cp_molar / M;  % J/(kg·K)
end

function cpR = polyCpR(T, a)
% Evaluate Cp/R for NASA7 (only first 5 coefficients used)
    cpR = a(1) + a(2).*T + a(3).*T.^2 + a(4).*T.^3 + a(5).*T.^4;
end
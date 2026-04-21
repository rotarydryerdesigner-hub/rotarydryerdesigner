function [Cp_gas, k_gas, rho_gas, mu_gas, Pr_gas, eps_gas] = ...
         gasPropsFromAirEquiv(Tg_C, P_Pa, fuelType)

    if nargin < 3
        fuelType = 'NaturalGas';
    end

    fuelType = strtrim(fuelType);

    % 1) Propiedades de aire "base" a Tg_C
    [k_air, rho_air, Cp_air, mu_air, nu_air, alpha_air, Pr_air, ~] = airProps(Tg_C);

    % 2) Factores de corrección según tipo de combustible
    switch fuelType
        case 'NaturalGas'
            fCp = 1.10;
            fk  = 0.95;
            fMu = 1.10;
            eps_gas = 0.90;

        case 'LPG'
            fCp = 1.12;
            fk  = 0.94;
            fMu = 1.12;
            eps_gas = 0.92;

        case 'Diesel'
            fCp = 1.14;
            fk  = 0.92;
            fMu = 1.18;
            eps_gas = 0.94;

        case 'FuelOil'
            fCp = 1.15;
            fk  = 0.90;
            fMu = 1.20;
            eps_gas = 0.95;

        case 'Biomass'
            fCp = 1.13;
            fk  = 0.93;
            fMu = 1.15;
            eps_gas = 0.94;

        case 'Coal'
            fCp = 1.16;
            fk  = 0.91;
            fMu = 1.22;
            eps_gas = 0.96;

        otherwise
            % genérico "flue gas"
            fCp = 1.10;
            fk  = 0.95;
            fMu = 1.10;
            eps_gas = 0.90;
    end

    % 3) Aplicar factores
    Cp_gas = fCp * Cp_air;   % J/kg·K
    k_gas  = fk  * k_air;    % W/m·K
    mu_gas = fMu * mu_air;   % Pa·s

    % 4) Densidad del gas (ideal, similar a aire)
    R_da    = 287.055;               
    T_K     = Tg_C + 273.15;
    rho_gas = P_Pa / (R_da * T_K);   

    % 5) Prandtl del gas
    Pr_gas = mu_gas * Cp_gas / k_gas;
end
function out = computeDesignGasFlow(in)

    % -------------------------
    % Required input fields
    % -------------------------
    req = { ...
        'Ms', 'X_in_db', 'X_target_db', ...
        'Tgin_C', 'Ts_in_C', 'Ts_out_C', ...
        'lossFrac', 'Cps', 'Yin', ...
        'Cpg', 'K_gas', 'rho_g', 'mu_g', ...
        'v_g_des'};

    for k = 1:numel(req)
        if ~isfield(in, req{k})
            error('Missing input field: %s', req{k});
        end
    end

    % -------------------------
    % Basic validation
    % -------------------------
    if ~isfinite(in.Ms) || in.Ms <= 0
        error('Wet feed mass flow must be positive.');
    end

    if ~isfinite(in.X_in_db) || in.X_in_db < 0
        error('Initial dry-basis moisture must be non-negative.');
    end

    if ~isfinite(in.X_target_db) || in.X_target_db < 0
        error('Target dry-basis moisture must be non-negative.');
    end

    if in.X_target_db >= in.X_in_db
        error(['Target dry-basis moisture must be lower than ', ...
               'initial dry-basis moisture.']);
    end

    temp_vec = [ ...
        in.Tgin_C, ...
        in.Ts_in_C, ...
        in.Ts_out_C];

    if any(~isfinite(temp_vec))
        error('Temperatures must be finite.');
    end

    if ~isfinite(in.Cps) || in.Cps <= 0
        error('Solid specific heat capacity Cps must be positive.');
    end

    if ~isfinite(in.Cpg) || in.Cpg <= 0
        error('Dry-gas specific heat capacity Cpg must be positive.');
    end

    gas_vec = [ ...
        in.K_gas, ...
        in.rho_g, ...
        in.mu_g];

    if any(~isfinite(gas_vec)) || any(gas_vec <= 0)
        error(['Gas properties are not defined correctly. ', ...
               'Check K_gas, rho_g, and mu_g.']);
    end

    if ~isfinite(in.Yin) || in.Yin < 0
        error('Inlet gas humidity ratio must be non-negative.');
    end

    if ~isfinite(in.v_g_des) || in.v_g_des <= 0
        error('Design gas velocity must be positive.');
    end

    if ~isfinite(in.lossFrac)
        error('Heat-loss fraction must be finite.');
    end

    % -------------------------
    % Heat-loss fraction
    % -------------------------
    lossFrac = in.lossFrac;

    % Allows user input as percentage, e.g. 5 instead of 0.05.
    if lossFrac > 1 && lossFrac <= 100
        lossFrac = lossFrac / 100;
    end

    if lossFrac < 0 || lossFrac >= 1
        error('Heat-loss fraction must be in the range [0, 1).');
    end

    % -------------------------
    % Solid-moisture convention
    % -------------------------
    % All solid moisture values use DRY BASIS:
    % X = kg water / kg dry solid.

    X_in_db     = in.X_in_db;
    X_target_db = in.X_target_db;

    % in.Ms is total wet-feed mass flow [kg/s].
    % For dry-basis moisture:
    %
    % Ms_wet = Ms_dry * (1 + X_in_db)
    %
    % Therefore:
    Ms_dry = in.Ms / (1 + X_in_db);

    % -------------------------
    % Thermophysical constants
    % -------------------------
    cp_s = in.Cps / 1000;     % kJ/(kg K)
    cp_g = in.Cpg / 1000;     % kJ/(kg K)

    if isfield(in, 'Cpv') && ...
       isfinite(in.Cpv) && ...
       in.Cpv > 0

        cp_v = in.Cpv / 1000; % kJ/(kg K)

    else

        cp_v = 1.884;         % kJ/(kg K), water vapor

    end

    cp_w    = 4.186;          % kJ/(kg K), liquid water
    lambda0 = 2502.3;         % kJ/kg, latent heat at 0 °C

    % -------------------------
    % Main mass calculations
    % -------------------------
    out.X_in_db     = X_in_db;
    out.X_target_db = X_target_db;

    % Keep these fields because the current App uses them.
    out.X1_db = X_in_db;
    out.X2_db = X_target_db;

    out.Ms_h  = in.Ms * 3600;      % kg/h wet feed
    out.Mds_h = Ms_dry * 3600;     % kg/h dry solid
    out.Ss_h  = out.Mds_h;         % kg/h dry solid

    out.m_evap_h = ...
        out.Ss_h * ...
        (X_in_db - X_target_db);    % kg/h evaporated water

    % -------------------------
    % Gas outlet temperature
    % -------------------------
    if isfield(in, 'Tg_out_C') && ...
       isfinite(in.Tg_out_C)

        Tg_out_C = in.Tg_out_C;

    else

        Tg_out_C = in.Ts_out_C;

        warning(['Tg_out_C was not provided. The code assumed ', ...
                 'Tg_out_C = Ts_out_C. This assumption should be justified.']);

    end

    if Tg_out_C >= in.Tgin_C
        error('Outlet gas temperature must be lower than inlet gas temperature.');
    end

    out.Tg_out_C_assumed = Tg_out_C;

    % -------------------------
    % Enthalpy calculations
    % -------------------------

    % Gas inlet enthalpy [kJ/kg dry gas]
    out.Hg_in = ...
        cp_g * in.Tgin_C + ...
        in.Yin * ...
        (lambda0 + cp_v * in.Tgin_C);

    % Gas outlet base enthalpy
    % without evaporated-water addition
    Hg_out_base = ...
        cp_g * Tg_out_C + ...
        in.Yin * ...
        (lambda0 + cp_v * Tg_out_C);

    out.Hg_out_base = Hg_out_base;

    % Solid-phase enthalpy [kJ/kg dry solid]
    Hs_in = ...
        cp_s * in.Ts_in_C + ...
        X_in_db * ...
        cp_w * in.Ts_in_C;

    Hs_out = ...
        cp_s * in.Ts_out_C + ...
        X_target_db * ...
        cp_w * in.Ts_out_C;

    out.Hs_in  = Hs_in;
    out.Hs_out = Hs_out;

    % Solid sensible and retained-moisture
    % sensible requirement [kJ/h]
    out.rhs_h = ...
        out.Ss_h * ...
        (Hs_out - Hs_in);

    % Energy required by evaporated water
    % leaving with outlet gas [kJ/h]
    out.evap_term_h = ...
        out.m_evap_h * ...
        (lambda0 + cp_v * Tg_out_C);

    % Total process heat load [kJ/h]
    out.Q_process_h = ...
        out.rhs_h + ...
        out.evap_term_h;

    if ~isfinite(out.Q_process_h) || ...
       out.Q_process_h <= 0

        error(['The calculated process heat load is non-positive or invalid. ', ...
               'Check temperatures, moisture contents, and heat capacities.']);

    end

    % Useful gas enthalpy drop [kJ/kg dry gas]
    out.deltaHg = ...
        out.Hg_in - ...
        Hg_out_base;

    if ~isfinite(out.deltaHg) || ...
       out.deltaHg <= 0

        error(['The gas enthalpy drop is non-positive. Check gas inlet ', ...
               'temperature, outlet gas temperature, and inlet humidity.']);

    end

    % -------------------------
    % Heat-loss model
    % -------------------------
    if isfield(in, 'lossBasis')

        lossBasis = ...
            lower(strtrim(in.lossBasis));

    else

        lossBasis = 'useful_duty';

    end

    switch lossBasis

        case {'useful_duty', 'useful', 'delta_h'}

            % Recommended formulation:
            % losses applied to useful gas enthalpy drop.
            Acoef = ...
                (1 - lossFrac) * ...
                out.deltaHg;

        case {'inlet_enthalpy', 'inlet'}

            % Legacy formulation.
            Acoef = ...
                (1 - lossFrac) * ...
                out.Hg_in - ...
                Hg_out_base;

        otherwise

            error(['Unknown lossBasis. Use ', ...
                   '''useful_duty'' or ''inlet_enthalpy''.']);

    end

    out.lossBasis = lossBasis;
    out.Acoef      = Acoef;

    if ~isfinite(Acoef) || Acoef <= 0

        error(['Acoef <= 0. The nominal balance does not close with the ', ...
               'current inputs. Check gas inlet/outlet temperatures, ', ...
               'loss fraction, and inlet humidity.']);

    end

    % -------------------------
    % Dry-gas mass flow
    % -------------------------
    out.GS_h = ...
        out.Q_process_h / ...
        Acoef;                      % kg dry gas/h

    if ~isfinite(out.GS_h) || ...
       out.GS_h <= 0

        error('Computed GS_h is non-positive or invalid.');

    end

    out.m_dot_g = ...
        out.GS_h / 3600;            % kg/s dry gas

    % -------------------------
    % Outlet gas humidity
    % -------------------------
    % Y = kg water/kg dry gas

    out.Yout = ...
        in.Yin + ...
        out.m_evap_h / out.GS_h;

    if ~isfinite(out.Yout) || ...
       out.Yout < in.Yin

        error('Computed outlet humidity ratio is invalid.');

    end

    % -------------------------
    % Saturation check
    % -------------------------
    if isfield(in, 'P_kPa') && ...
       isfinite(in.P_kPa) && ...
       in.P_kPa > 0

        P_kPa = in.P_kPa;

    else

        P_kPa = 101.325;

    end

    out.P_kPa    = P_kPa;
    out.Ysat_out = NaN;

    % Saturation check is meaningful for typical
    % outlet temperatures below the boiling point
    % at near-atmospheric pressure.
    if Tg_out_C > -40 && Tg_out_C < 100

        Pws_kPa = ...
            localWaterSaturationPressure_kPa(Tg_out_C);

        if Pws_kPa < P_kPa

            out.Ysat_out = ...
                0.62198 * ...
                Pws_kPa / ...
                (P_kPa - Pws_kPa);

            if out.Yout > out.Ysat_out

                error(['The computed outlet gas is supersaturated: ', ...
                       'Yout > Ysat(Tg_out_C). Increase dry-gas flow, ', ...
                       'increase outlet gas temperature, or revise ', ...
                       'the target moisture.']);

            end
        end
    end

    % -------------------------
    % Preliminary gas-flow geometry
    % -------------------------
    out.v_g_des = in.v_g_des;

    out.rho_g = in.rho_g;
    out.mu_g  = in.mu_g;
    out.K_gas = in.K_gas;
    out.Cpg   = in.Cpg;

    out.A_gas_m2 = ...
        out.m_dot_g / ...
        (in.rho_g * in.v_g_des);

    out.D_eq_m = ...
        sqrt(4 * out.A_gas_m2 / pi);

    out.Re_g = ...
        in.rho_g * ...
        in.v_g_des * ...
        out.D_eq_m / ...
        in.mu_g;

    out.lossFrac = lossFrac;

    % Solid-moisture convention used by this function.
    out.moistureBasis = 'db';

end


function Pws_kPa = localWaterSaturationPressure_kPa(T_C)
% localWaterSaturationPressure_kPa
% Saturation pressure of water vapor using a Buck-type equation.
% Valid with acceptable engineering accuracy for common drying outlet
% temperatures below 100 °C.

    Pws_kPa = ...
        0.61121 * ...
        exp((18.678 - T_C/234.5) * ...
        (T_C / (257.14 + T_C)));

end
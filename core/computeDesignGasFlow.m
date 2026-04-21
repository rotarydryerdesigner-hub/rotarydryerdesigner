function out = computeDesignGasFlow(in)
% computeDesignGasFlow
% Calculates the nominal dry-gas flow required from a simplified mass and
% energy balance for rotary dryer preliminary design.
%
% Input:
%   in : struct with required fields
%
% Output:
%   out : struct with calculated variables

    % -------------------------
    % Validate required inputs
    % -------------------------
    req = { ...
        'Ms', 'hum', 'h_target', 'Tgin_C', 'Ts_in_C', 'Ts_out_C', ...
        'lossFrac', 'Cps', 'Yin', 'Cpg', 'K_gas', 'rho_g', 'mu_g', ...
        'v_g_des'};

    for k = 1:numel(req)
        if ~isfield(in, req{k})
            error('Missing input field: %s', req{k});
        end
    end

    gas_vec = [in.Cpg, in.K_gas, in.rho_g, in.mu_g, in.Yin];
    if any(~isfinite(gas_vec)) || any(gas_vec <= 0)
        error(['Gas properties are not defined correctly. ', ...
               'Check Cpg, K_gas, rho_g, mu_g, and Yin.']);
    end

    if ~isfinite(in.v_g_des) || in.v_g_des <= 0
        error('Design gas velocity must be positive.');
    end

    % -------------------------
    % Main calculations
    % -------------------------
    Cpv = Cp_h2o(in.Ts_out_C);   % J/(kg·K)

    out.X1_db = in.hum;
    out.X2_db = in.h_target;

    out.Ms_h = in.Ms * 3600;     % kg/h wet feed
    out.Ss_h = out.Ms_h;         % mantiene tu misma suposición actual
    out.m_evap_h = out.Ss_h * (out.X1_db - out.X2_db);

    cp_s = in.Cps / 1000;        % kJ/kg/K
    cp_w = Cpv / 1000;           % kJ/kg/K

    out.Hg_in = (1.005 + 1.884*in.Yin)*in.Tgin_C + 2502.3*in.Yin;

    Hs_in     = cp_s*in.Ts_in_C  + out.X1_db*cp_w*in.Ts_in_C;
    out.Hs_out = cp_s*in.Ts_out_C + out.X2_db*cp_w*in.Ts_out_C;

    out.rhs_h = out.Ss_h * (out.Hs_out - Hs_in);

    Tg_out_C = in.Ts_out_C;
    Hg_out_sensible = (1.005 + 1.884*in.Yin)*Tg_out_C;

    Acoef = out.Hg_in*(1 - in.lossFrac) - Hg_out_sensible - 2502.3*in.Yin;

    if ~isfinite(Acoef) || Acoef <= 0
        error(['Acoef <= 0. The nominal balance does not close with the ', ...
               'current inputs.']);
    end

    out.GS_h = (out.rhs_h + 2502.3*out.m_evap_h) / Acoef;

    if ~isfinite(out.GS_h) || out.GS_h <= 0
        error('Computed GS_h is non-positive or invalid.');
    end

    out.Yout = in.Yin + out.m_evap_h / out.GS_h;
    out.m_dot_g = out.GS_h / 3600;
    out.v_g_des = in.v_g_des;
end
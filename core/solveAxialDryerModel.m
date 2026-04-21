function out = solveAxialDryerModel(in)
% solveAxialDryerModel
% Solves the one-dimensional axial rotary dryer model and returns profiles
% and summary quantities needed by the App.

    % -------------------------
    % Required input fields
    % -------------------------
    req = { ...
        'Ms', 'hum', 'h_target', ...
        'Tgin_C', 'Ts_in_C', 'Ts_out_C', ...
        'lossFrac', 'f_fill', 'Vp', 'phiD0', ...
        'Tamb_C', 'h_alt', ...
        'den_s', 'Ks', 'Cps', 'es', ...
        'h_step', 'Lmax', 'H_tol', ...
        'A', 'E', 'Heq', ...
        'Yin', 'Cpg', 'K_gas', 'rho_g', 'mu_g', 'eg', ...
        'Di', 'Mgas', ...
        'Ss_h', 'Ms_h', 'm_evap_h', 'rhs_h'};

    for k = 1:numel(req)
        if ~isfield(in, req{k})
            error('Missing input field: %s', req{k});
        end
    end

    % -------------------------
    % Basic validation
    % -------------------------
    gas_vec = [in.Cpg, in.K_gas, in.rho_g, in.mu_g, in.Yin];
    if any(~isfinite(gas_vec)) || any(gas_vec <= 0)
        error(['Gas properties not defined or non-positive. ', ...
               'Check Cpg, K_gas, rho_g, mu_g, and Yin.']);
    end

    if ~isfinite(in.Di) || in.Di <= 0
        error('Dryer diameter must be positive.');
    end

    if ~isfinite(in.Mgas) || in.Mgas <= 0
        error('Gas mass flow must be positive.');
    end

    if ~isfinite(in.h_step) || in.h_step <= 0
        error('Axial step must be positive.');
    end

    if ~isfinite(in.Lmax) || in.Lmax <= 0
        error('Maximum axial length must be positive.');
    end

    % -------------------------
    % Constants
    % -------------------------
    CBosm = 5.670374419e-8;
    R     = 8.314;
    Hv    = 2.22e6;
    ew    = 0.9;
    ewo   = 0.9;

      % -------------------------
    % Derived gas properties
    % -------------------------
    nu_g    = in.mu_g / in.rho_g;
    alpha_g = in.K_gas / (in.rho_g * in.Cpg);
    Pr_g    = nu_g / alpha_g;

    Cpv = Cp_h2o(in.Ts_out_C);

    % -------------------------
    % Geometry and kinematics
    % -------------------------
    Ri   = in.Di / 2;
    Asec = pi * (in.Di^2) / 4;

    v_gas = in.Mgas / (in.rho_g * Asec);
    RPM   = 60 * in.Vp / (pi * in.Di);

    Ss      = in.Ss_h / 3600;
    rho_s   = in.den_s;
    N_rev_s = RPM / 60;

    phi = in.phiD0;
    if phi > 1
        phi = phi / 100;
    end

    S_design  = 0.3344 * Ss / (max(phi, eps) * rho_s * (max(N_rev_s, eps)^0.9) * max(in.Di, eps));
    Ss_design = 100 * S_design;

    Vs = in.Ms / max(in.den_s * in.f_fill * Asec, eps);

    betta = ((3 * in.Ms) / (Vs * in.den_s * (Ri^2)))^(1/3);

    Lwcs  = 2 * Ri * betta;
    Lwncs = Ri * (2*pi - 2*betta);
    Lss   = 2 * Ri * sin(betta);

    De   = (0.5 * in.Di * ((2*pi) - 2*betta + sin(2*betta))) / (pi - betta + sin(betta));
    Asnc = (pi * (De^2)) / 4;

    % -------------------------
    % Convective coefficients
    % -------------------------
    Re_gas = (v_gas * in.Di) / nu_g;

    if Re_gas >= 0.4 && Re_gas < 4
        Cc = 0.989; mc = 0.33;
    elseif Re_gas >= 4 && Re_gas < 40
        Cc = 0.911; mc = 0.385;
    elseif Re_gas >= 40 && Re_gas < 4000
        Cc = 0.686; mc = 0.466;
    elseif Re_gas >= 4000 && Re_gas < 40000
        Cc = 0.196; mc = 0.618;
    else
        Cc = 0.027; mc = 0.805;
    end

    Nu_cwo = Cc * (Re_gas^mc) * Pr_g^(1/3);
    ho     = (Nu_cwo * in.K_gas) / in.Di;

    hgw = Cc * Pr_g^(1/3) * Re_gas^mc;
    hgs = hgw * 0.10;

    Ab  = in.Ks / (in.den_s * in.Cps);
    hsw = 0.116 * in.Ks * (((RPM/60) * (Ri^2) * 2 * betta) / Ab)^0.3 / in.Di;

    % -------------------------
    % Initial conditions
    % -------------------------
    Tair_K = in.Tamb_C + 273.15;
    Ts0    = in.Ts_in_C + 273.15;
    Tg0    = in.Tgin_C  + 273.15;

    K1_0   = (hgs + hgw) * Asnc;
    K3_0   = hsw * Asnc;
    hcwo_0 = ho * Asnc;

    Tw0 = (K1_0*Tg0 + K3_0*Ts0 + hcwo_0*Tair_K) / (K1_0 + K3_0 + hcwo_0);
    Two_K = Tw0;

    maxN = ceil(in.Lmax / in.h_step) + 1;

    x_ax  = zeros(1, maxN);
    Tg_ax = zeros(1, maxN);
    Ts_ax = zeros(1, maxN);
    Tw_ax = zeros(1, maxN);
    H_ax  = zeros(1, maxN);
    Y_ax  = zeros(1, maxN);

    i = 1;
    x_ax(i)  = 0;
    Tg_ax(i) = Tg0;
    Ts_ax(i) = Ts0;
    Tw_ax(i) = Tw0;
    H_ax(i)  = in.hum;
    Y_ax(i)  = in.Yin;

    % -------------------------
    % Kinetics and ODEs
    % -------------------------
    Kc_fun = @(Tg) in.A * exp(-in.E ./ (R * Tg));

    Ms_per_m = in.Ms / Vs;

    F_xH = @(Hloc, TgLoc) -(Kc_fun(TgLoc) .* (Hloc - in.Heq)) / Vs;

    F_xY = @(Yloc, Hloc, TgLoc) ...
        (Kc_fun(TgLoc) .* (Hloc - in.Heq)) * (in.Ms / max(in.Mgas, eps)) / v_gas;

    F_xTg = @(TgLoc, TsLoc, TwLoc, K1loc, K2loc, Raloc) ...
        ((-K1loc*(TgLoc - TsLoc) - K2loc*(TgLoc - TwLoc) + ...
        Raloc*Ms_per_m*Cpv*(TsLoc - TgLoc)) / (max(in.Mgas, eps) * in.Cpg));

    F_xTs = @(TgLoc, TsLoc, TwLoc, K1loc, K3loc, Raloc) ...
        ((K1loc*(TgLoc - TsLoc) + K3loc*(TwLoc - TsLoc) - Hv*Raloc*in.Ms) / ...
        (in.Cps * max(in.Ms, eps)));

    stopIdx = maxN;


    while i < maxN

        hwa  = ewo * CBosm * (Two_K^3) * ...
            (1 + (Tair_K/Two_K) + (Tair_K/Two_K)^2 + (Tair_K/Two_K)^3);
        hcwo = ho + hwa;

        hgw_r = in.eg * ew * CBosm * (Tg_ax(i)^3) * ...
            (1 + (Tw_ax(i)/Tg_ax(i)) + (Tw_ax(i)/Tg_ax(i))^2 + (Tw_ax(i)/Tg_ax(i))^3);

        hgs_r = in.es * in.eg * CBosm * (Tg_ax(i)^3) * ...
            (1 + (Ts_ax(i)/Tg_ax(i)) + (Ts_ax(i)/Tg_ax(i))^2 + (Ts_ax(i)/Tg_ax(i))^3);

        dPM = (Tw_ax(i) - Ts_ax(i));
        if abs(dPM) < 1e-9
            dPM = sign(dPM + 1e-12) * 1e-9;
        end

        hws_r = (CBosm * in.es * ew * ((Tw_ax(i)^4) - (Ts_ax(i)^4))) / dPM;

        K1 = (hgs + hgs_r) * Lss;
        K2 = (hgw + hgw_r) * Lwncs;
        K3 = (hsw * Lwcs + hws_r * Lss);

        Ra = Kc_fun(Tg_ax(i)) * (H_ax(i) - in.Heq);
        Ra = max(Ra, 0);

        Qin_s = K1 * (Tg_ax(i) - Ts_ax(i)) + K3 * (Tw_ax(i) - Ts_ax(i));
        Qin_s = max(Qin_s, 0);
        Ra_max = Qin_s / max(Hv * Ms_per_m, eps);
        Ra = min(Ra, Ra_max);

        H1  = in.h_step * F_xH(H_ax(i), Tg_ax(i));
        Y1  = in.h_step * F_xY(Y_ax(i), H_ax(i), Tg_ax(i));
        Tg1 = in.h_step * F_xTg(Tg_ax(i), Ts_ax(i), Tw_ax(i), K1, K2, Ra);
        Ts1 = in.h_step * F_xTs(Tg_ax(i), Ts_ax(i), Tw_ax(i), K1, K3, Ra);

        H2  = H_ax(i)  + 0.5 * H1;
        Y2l = Y_ax(i)  + 0.5 * Y1;
        Tg2 = Tg_ax(i) + 0.5 * Tg1;
        Ts2 = Ts_ax(i) + 0.5 * Ts1;
        Tw2 = (K1*Tg2 + K3*Ts2 + hcwo*Tair_K) / (K1 + K3 + hcwo);

        H2k  = in.h_step * F_xH(H2, Tg2);
        Y2k  = in.h_step * F_xY(Y2l, H2, Tg2);
        Tg2k = in.h_step * F_xTg(Tg2, Ts2, Tw2, K1, K2, Ra);
        Ts2k = in.h_step * F_xTs(Tg2, Ts2, Tw2, K1, K3, Ra);

        H3  = H_ax(i)  + 0.5 * H2k;
        Y3l = Y_ax(i)  + 0.5 * Y2k;
        Tg3 = Tg_ax(i) + 0.5 * Tg2k;
        Ts3 = Ts_ax(i) + 0.5 * Ts2k;
        Tw3 = (K1*Tg3 + K3*Ts3 + hcwo*Tair_K) / (K1 + K3 + hcwo);

        H3k  = in.h_step * F_xH(H3, Tg3);
        Y3k  = in.h_step * F_xY(Y3l, H3, Tg3);
        Tg3k = in.h_step * F_xTg(Tg3, Ts3, Tw3, K1, K2, Ra);
        Ts3k = in.h_step * F_xTs(Tg3, Ts3, Tw3, K1, K3, Ra);

        H4  = H_ax(i)  + H3k;
        Y4l = Y_ax(i)  + Y3k;
        Tg4 = Tg_ax(i) + Tg3k;
        Ts4 = Ts_ax(i) + Ts3k;
        Tw4 = (K1*Tg4 + K3*Ts4 + hcwo*Tair_K) / (K1 + K3 + hcwo);

        H4k  = in.h_step * F_xH(H4, Tg4);
        Y4k  = in.h_step * F_xY(Y4l, H4, Tg4);
        Tg4k = in.h_step * F_xTg(Tg4, Ts4, Tw4, K1, K2, Ra);
        Ts4k = in.h_step * F_xTs(Tg4, Ts4, Tw4, K1, K3, Ra);

        H_next  = H_ax(i)  + (1/6) * (H1  + 2*H2k  + 2*H3k  + H4k);
        Y_next  = Y_ax(i)  + (1/6) * (Y1  + 2*Y2k  + 2*Y3k  + Y4k);
        Tg_next = Tg_ax(i) + (1/6) * (Tg1 + 2*Tg2k + 2*Tg3k + Tg4k);
        Ts_next = Ts_ax(i) + (1/6) * (Ts1 + 2*Ts2k + 2*Ts3k + Ts4k);

        Tw_next = (K1*Tg_next + K3*Ts_next + hcwo*Tair_K) / (K1 + K3 + hcwo);

        if (H_ax(i) > in.h_target) && ...
           (H_next <= in.h_target || abs(H_next - in.h_target) <= in.H_tol)

            alpha = 1;
            if H_next ~= H_ax(i)
                alpha = (in.h_target - H_ax(i)) / (H_next - H_ax(i));
                alpha = max(0, min(1, alpha));
            end

            i = i + 1;
            x_ax(i)  = x_ax(i-1) + alpha * in.h_step;

            H_ax(i)  = in.h_target;
            Y_ax(i)  = Y_ax(i-1)  + alpha * (Y_next  - Y_ax(i-1));
            Tg_ax(i) = Tg_ax(i-1) + alpha * (Tg_next - Tg_ax(i-1));
            Ts_ax(i) = Ts_ax(i-1) + alpha * (Ts_next - Ts_ax(i-1));
            Tw_ax(i) = Tw_ax(i-1) + alpha * (Tw_next - Tw_ax(i-1));

            stopIdx = i;
            break
        end

        i = i + 1;
        x_ax(i)  = x_ax(i-1) + in.h_step;

        H_ax(i)  = H_next;
        Y_ax(i)  = Y_next;
        Tg_ax(i) = Tg_next;
        Ts_ax(i) = Ts_next;
        Tw_ax(i) = Tw_next;
    end

    % -------------------------
    % Post-processing
    % -------------------------
    x_ax  = x_ax(1:stopIdx);
    Tg_ax = Tg_ax(1:stopIdx);
    Ts_ax = Ts_ax(1:stopIdx);
    H_ax  = H_ax(1:stopIdx);
    Tw_ax = Tw_ax(1:stopIdx);
    Y_ax  = Y_ax(1:stopIdx);

    L_req = x_ax(end);

    Tg_out_K       = Tg_ax(end);
    Ts_out_K_model = Ts_ax(end);

    Ac     = pi * in.Di^2 / 4;
    Ms_dry = in.Ms * (1 - in.hum);

    Trx = (L_req * phi * in.den_s * Ac) / max(Ms_dry, eps) / 60;
    Trx_min = Trx / 60;

    % -------------------------
    % Output structure
    % -------------------------
    out = struct();

    out.x_ax = x_ax;
    out.Tg_ax = Tg_ax;
    out.Ts_ax = Ts_ax;
    out.Tw_ax = Tw_ax;
    out.H_ax = H_ax;
    out.Y_ax = Y_ax;

    out.L_req = L_req;
    out.Tg_out_K = Tg_out_K;
    out.Ts_out_K_model = Ts_out_K_model;

    out.RPM = RPM;
    out.Ss_design = Ss_design;
    out.Trx = Trx;
    out.Trx_min = Trx_min;

    out.phi = phi;
    out.v_gas = v_gas;
    out.Re_gas = Re_gas;

    out.Cpv = Cpv;
    out.Pr_g = Pr_g;
    out.nu_g = nu_g;
    out.alpha_g = alpha_g;

    out.ho = ho;
    out.hgw = hgw;
    out.hgs = hgs;
    out.hsw = hsw;

    out.solverOK = true;
    out.targetOK = (H_ax(end) <= in.h_target + in.H_tol);
    out.validityWarning = (Ss_design > 10);
end


    


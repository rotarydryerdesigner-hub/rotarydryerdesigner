function out = solveAxialDryerModel(in)
% solveAxialDryerModel
% Solves the one-dimensional axial rotary dryer model and returns profiles
% and summary quantities needed by the App.
% Solid moisture is handled internally on a dry basis:
% X = kg water / kg dry solid.

    % -------------------------
    % Required input fields
    % -------------------------
    req = { ...
        'Ms', 'X_in_db', 'X_target_db', ...
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
    % Solid-moisture convention
    % -------------------------
    % All solid-moisture calculations in this solver use DRY BASIS:
    % X = kg water / kg dry solid.
    X_in_db     = in.X_in_db;
    X_target_db = in.X_target_db;
    Xeq_db      = in.Heq;
    X_tol       = in.H_tol;

    if any(~isfinite([X_in_db, X_target_db, Xeq_db, X_tol]))
        error('Solid-moisture inputs must be finite.');
    end

    if X_in_db < 0 || X_target_db < 0 || Xeq_db < 0 || X_tol < 0
        error('Solid-moisture values and tolerance must be non-negative.');
    end

    if X_target_db >= X_in_db
        error('Target dry-basis moisture must be lower than inlet dry-basis moisture.');
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

    S_design = 0.3344 * Ss / ...
        (max(phi, eps) * rho_s * ...
        (max(N_rev_s, eps)^0.9) * max(in.Di, eps));

    Ss_design = 100 * S_design;

    Vs = in.Ms / ...
        max(in.den_s * in.f_fill * Asec, eps);

    betta = ((3 * in.Ms) / ...
        (Vs * in.den_s * (Ri^2)))^(1/3);

    Lwcs  = 2 * Ri * betta;
    Lwncs = Ri * (2*pi - 2*betta);
    Lss   = 2 * Ri * sin(betta);

    De = (0.5 * in.Di * ...
        ((2*pi) - 2*betta + sin(2*betta))) / ...
        (pi - betta + sin(betta));

    Asnc = (pi * (De^2)) / 4;

    % -------------------------
    % Convective coefficients
    % -------------------------
    Re_gas = (v_gas * in.Di) / nu_g;

    Cx=0.4;
    Bx=0.54;


    if Re_gas >= 0.4 && Re_gas < 4
        Cc = 0.989;
        mc = 0.33;

    elseif Re_gas >= 4 && Re_gas < 40
        Cc = 0.911;
        mc = 0.385;

    elseif Re_gas >= 40 && Re_gas < 4000
        Cc = 0.686;
        mc = 0.466;

    elseif Re_gas >= 4000 && Re_gas < 40000
        Cc = 0.196;
        mc = 0.618;

    else
        Cc = 0.027;
        mc = 0.805;
    end

    Nu_cwo = Cc * (Re_gas^mc) * Pr_g^(1/3);

    ho = (Nu_cwo * in.K_gas) / in.Di;


 

% -------------------------
% Gas-solid coefficient
% Tscheng and Watkinson correlation
% -------------------------

G_g = (in.Mgas*3600)/ Asnc;       % [kg/(m^2 h)]


hgs =  Cx* G_g^Bx;                % [W/(m^2 K)]
      

Nu_gw = Cc * (Re_gas^mc) * Pr_g^(1/3);

%hgw = (Nu_gw * in.K_gas) / in.Di;   % [W/(m^2 K)]


hgw =  0.098* G_g^0.89;                % [W/(m^2 K)]


Ab = in.Ks / (in.den_s * in.Cps);

    hsw = 0.116 * in.Ks * ...
        (((RPM/60) * (Ri^2) * 2 * betta) / Ab)^0.3 / in.Di;

    % -------------------------
    % Initial conditions
    % -------------------------
    Tair_K = in.Tamb_C + 273.15;
    Ts0    = in.Ts_in_C + 273.15;
    Tg0    = in.Tgin_C + 273.15;
    Two_0=50+273;
    
     
      hwa = ewo * CBosm * (Two_0^3) * ...
            (1 + ...
            (Tair_K/Two_0) + ...
            (Tair_K/Two_0)^2 + ...
            (Tair_K/Two_0)^3);


    K1_0   = (hgs + hgw) * Asnc;
    K3_0   = hsw * Asnc;
    hcwo_0 = (ho+ hwa) * pi * in.Di;




    
    Tw0 = ...
        (K1_0*Tg0 + ...
         K3_0*Ts0 + ...
         hcwo_0*Tair_K) / ...
        (K1_0 + K3_0 + hcwo_0);

    Two_K = Tw0;

    maxN = ceil(in.Lmax / in.h_step) + 1;

    x_ax  = zeros(1,maxN);
    Tg_ax = zeros(1,maxN);
    Ts_ax = zeros(1,maxN);
    Tw_ax = zeros(1,maxN);
    X_ax  = zeros(1,maxN);
    Y_ax  = zeros(1,maxN);

    i = 1;

    x_ax(i)  = 0;
    Tg_ax(i) = Tg0;
    Ts_ax(i) = Ts0;
    Tw_ax(i) = Tw0;

    % Solid moisture: dry basis
    X_ax(i) = X_in_db;

    Y_ax(i) = in.Yin;

    % -------------------------
    % Kinetics and ODEs
    % -------------------------
    Kc_fun = @(Tg) ...
        in.A * exp(-in.E ./ (R*Tg));

    Ms_per_m = in.Ms / Vs;

    % Solid moisture balance
    F_xX = @(Xloc,TgLoc) ...
        -(Kc_fun(TgLoc) .* ...
        (Xloc - Xeq_db)) / Vs;

    % Gas humidity balance
    F_xY = @(Yloc,Xloc,TgLoc) ...
        (Kc_fun(TgLoc) .* ...
        (Xloc - Xeq_db)) * ...
        (in.Ms / max(in.Mgas,eps)) / v_gas;

    % Gas temperature balance
    F_xTg = @(TgLoc,TsLoc,TwLoc,K1loc,K2loc,Raloc) ...
        ((-K1loc*(TgLoc-TsLoc) ...
        - K2loc*(TgLoc-TwLoc) ...
        + Raloc*Ms_per_m*Cpv*(TsLoc-TgLoc)) / ...
        (max(in.Mgas,eps)*in.Cpg));

    % Solid temperature balance
    F_xTs = @(TgLoc,TsLoc,TwLoc,K1loc,K3loc,Raloc) ...
        ((K1loc*(TgLoc-TsLoc) ...
        + K3loc*(TwLoc-TsLoc) ...
        - Hv*Raloc*in.Ms) / ...
        (in.Cps*max(in.Ms,eps)));

    stopIdx = maxN;

    % -------------------------
    % Axial integration
    % -------------------------
    while i < maxN

        hwa = ewo * CBosm * (Two_K^3) * ...
            (1 + ...
            (Tair_K/Two_K) + ...
            (Tair_K/Two_K)^2 + ...
            (Tair_K/Two_K)^3);

        hcwo = ho + hwa;

        hgw_r = in.eg * ew * CBosm * ...
            (Tg_ax(i)^3) * ...
            (1 + ...
            (Tw_ax(i)/Tg_ax(i)) + ...
            (Tw_ax(i)/Tg_ax(i))^2 + ...
            (Tw_ax(i)/Tg_ax(i))^3);

        hgs_r = in.es * in.eg * CBosm * ...
            (Tg_ax(i)^3) * ...
            (1 + ...
            (Ts_ax(i)/Tg_ax(i)) + ...
            (Ts_ax(i)/Tg_ax(i))^2 + ...
            (Ts_ax(i)/Tg_ax(i))^3);

        dPM = Tw_ax(i) - Ts_ax(i);

        if abs(dPM) < 1e-9
            dPM = sign(dPM + 1e-12) * 1e-9;
        end

        hws_r = ...
            (CBosm * in.es * ew * ...
            ((Tw_ax(i)^4) - (Ts_ax(i)^4))) / dPM;

        % -------------------------
        % Effective heat-transfer terms
        % -------------------------
        K1 = (hgs + hgs_r) * Lss;

        K2 = (hgw + hgw_r) * Lwncs;

        K3 = ...
            (hsw * Lwcs + ...
            hws_r * Lss);

        % -------------------------
        % Drying rate
        % -------------------------
        Ra = ...
            Kc_fun(Tg_ax(i)) * ...
            (X_ax(i) - Xeq_db);

        Ra = max(Ra,0);

        Qin_s = ...
            K1*(Tg_ax(i)-Ts_ax(i)) + ...
            K3*(Tw_ax(i)-Ts_ax(i));

        Qin_s = max(Qin_s,0);

        Ra_max = ...
            Qin_s / max(Hv*Ms_per_m,eps);

        Ra = min(Ra,Ra_max);

        % =====================================================
        % RK4 - Stage 1
        % =====================================================
        X1 = ...
            in.h_step * ...
            F_xX(X_ax(i),Tg_ax(i));

        Y1 = ...
            in.h_step * ...
            F_xY(Y_ax(i),X_ax(i),Tg_ax(i));

        Tg1 = ...
            in.h_step * ...
            F_xTg( ...
            Tg_ax(i), ...
            Ts_ax(i), ...
            Tw_ax(i), ...
            K1,K2,Ra);

        Ts1 = ...
            in.h_step * ...
            F_xTs( ...
            Tg_ax(i), ...
            Ts_ax(i), ...
            Tw_ax(i), ...
            K1,K3,Ra);

        % =====================================================
        % RK4 - Stage 2
        % =====================================================
        X2  = X_ax(i)  + 0.5*X1;
        Y2l = Y_ax(i)  + 0.5*Y1;
        Tg2 = Tg_ax(i) + 0.5*Tg1;
        Ts2 = Ts_ax(i) + 0.5*Ts1;

        Tw2 = ...
            (K1*Tg2 + ...
            K3*Ts2 + ...
            hcwo*Tair_K) / ...
            (K1 + K3 + hcwo);

        X2k = ...
            in.h_step * ...
            F_xX(X2,Tg2);

        Y2k = ...
            in.h_step * ...
            F_xY(Y2l,X2,Tg2);

        Tg2k = ...
            in.h_step * ...
            F_xTg( ...
            Tg2,Ts2,Tw2, ...
            K1,K2,Ra);

        Ts2k = ...
            in.h_step * ...
            F_xTs( ...
            Tg2,Ts2,Tw2, ...
            K1,K3,Ra);

        % =====================================================
        % RK4 - Stage 3
        % =====================================================
        X3  = X_ax(i)  + 0.5*X2k;
        Y3l = Y_ax(i)  + 0.5*Y2k;
        Tg3 = Tg_ax(i) + 0.5*Tg2k;
        Ts3 = Ts_ax(i) + 0.5*Ts2k;

        Tw3 = ...
            (K1*Tg3 + ...
            K3*Ts3 + ...
            hcwo*Tair_K) / ...
            (K1 + K3 + hcwo);

        X3k = ...
            in.h_step * ...
            F_xX(X3,Tg3);

        Y3k = ...
            in.h_step * ...
            F_xY(Y3l,X3,Tg3);

        Tg3k = ...
            in.h_step * ...
            F_xTg( ...
            Tg3,Ts3,Tw3, ...
            K1,K2,Ra);

        Ts3k = ...
            in.h_step * ...
            F_xTs( ...
            Tg3,Ts3,Tw3, ...
            K1,K3,Ra);

        % =====================================================
        % RK4 - Stage 4
        % =====================================================
        X4  = X_ax(i)  + X3k;
        Y4l = Y_ax(i)  + Y3k;
        Tg4 = Tg_ax(i) + Tg3k;
        Ts4 = Ts_ax(i) + Ts3k;

        Tw4 = ...
            (K1*Tg4 + ...
            K3*Ts4 + ...
            hcwo*Tair_K) / ...
            (K1 + K3 + hcwo);

        X4k = ...
            in.h_step * ...
            F_xX(X4,Tg4);

        Y4k = ...
            in.h_step * ...
            F_xY(Y4l,X4,Tg4);

        Tg4k = ...
            in.h_step * ...
            F_xTg( ...
            Tg4,Ts4,Tw4, ...
            K1,K2,Ra);

        Ts4k = ...
            in.h_step * ...
            F_xTs( ...
            Tg4,Ts4,Tw4, ...
            K1,K3,Ra);

        % =====================================================
        % RK4 update
        % =====================================================
        X_next = ...
            X_ax(i) + ...
            (1/6)*( ...
            X1 + ...
            2*X2k + ...
            2*X3k + ...
            X4k);

        Y_next = ...
            Y_ax(i) + ...
            (1/6)*( ...
            Y1 + ...
            2*Y2k + ...
            2*Y3k + ...
            Y4k);

        Tg_next = ...
            Tg_ax(i) + ...
            (1/6)*( ...
            Tg1 + ...
            2*Tg2k + ...
            2*Tg3k + ...
            Tg4k);

        Ts_next = ...
            Ts_ax(i) + ...
            (1/6)*( ...
            Ts1 + ...
            2*Ts2k + ...
            2*Ts3k + ...
            Ts4k);

        Tw_next = ...
            (K1*Tg_next + ...
            K3*Ts_next + ...
            hcwo*Tair_K) / ...
            (K1 + K3 + hcwo);

        % -------------------------
        % Target-moisture detection
        % -------------------------
        if (X_ax(i) > X_target_db) && ...
           (X_next <= X_target_db || ...
           abs(X_next-X_target_db) <= X_tol)

            alpha = 1;

            if X_next ~= X_ax(i)

                alpha = ...
                    (X_target_db-X_ax(i)) / ...
                    (X_next-X_ax(i));

                alpha = max(0,min(1,alpha));

            end

            i = i + 1;

            x_ax(i) = ...
                x_ax(i-1) + ...
                alpha*in.h_step;

            X_ax(i) = X_target_db;

            Y_ax(i) = ...
                Y_ax(i-1) + ...
                alpha*(Y_next-Y_ax(i-1));

            Tg_ax(i) = ...
                Tg_ax(i-1) + ...
                alpha*(Tg_next-Tg_ax(i-1));

            Ts_ax(i) = ...
                Ts_ax(i-1) + ...
                alpha*(Ts_next-Ts_ax(i-1));

            Tw_ax(i) = ...
                Tw_ax(i-1) + ...
                alpha*(Tw_next-Tw_ax(i-1));

            stopIdx = i;

            break
        end

        i = i + 1;

        x_ax(i) = ...
            x_ax(i-1) + in.h_step;

        X_ax(i)  = X_next;
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
    X_ax  = X_ax(1:stopIdx);
    Tw_ax = Tw_ax(1:stopIdx);
    Y_ax  = Y_ax(1:stopIdx);

    L_req = x_ax(end);

    Tg_out_K = Tg_ax(end);

    Ts_out_K_model = Ts_ax(end);

    Ac = pi * in.Di^2 / 4;

    % Dry-solid mass flow
    Ms_dry = in.Ss_h / 3600;

    Trx = ...
        (L_req * phi * in.den_s * Ac) / ...
        max(Ms_dry,eps) / 60;

    Trx_min = Trx / 60;

    % -------------------------
    % Output structure
    % -------------------------
    out = struct();

    out.x_ax  = x_ax;
    out.Tg_ax = Tg_ax;
    out.Ts_ax = Ts_ax;
    out.Tw_ax = Tw_ax;

    % Canonical solid-moisture variable
    out.X_ax = X_ax;

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

    % -------------------------
    % Moisture-basis information
    % -------------------------
    out.X_in_db = X_in_db;

    out.X_target_db = X_target_db;

    out.Xeq_db = Xeq_db;

    out.moistureBasis = 'db';

    % -------------------------
    % Solver status
    % -------------------------
    out.solverOK = true;

    out.targetOK = ...
        (X_ax(end) <= X_target_db + X_tol);

    out.validityWarning = ...
        (Ss_design > 40);

end


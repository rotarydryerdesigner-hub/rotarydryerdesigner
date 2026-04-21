function fuel = compute_fuel_consumption(fuelType, rhs_h, m_evap_h, lossFrac, Ms_h, Ss_h)

switch lower(strtrim(fuelType))
    case 'natural gas'
        LHV_kJkg = 50000;
        eta_burner = 0.85;

    case 'lpg'
        LHV_kJkg = 46000;
        eta_burner = 0.85;

    case 'diesel'
        LHV_kJkg = 42800;
        eta_burner = 0.82;

    case 'fuel oil'
        LHV_kJkg = 40500;
        eta_burner = 0.80;

    case 'biomass pellets'
        LHV_kJkg = 17000;
        eta_burner = 0.75;

    case 'coal'
        LHV_kJkg = 25000;
        eta_burner = 0.78;

    otherwise
        error('Unsupported fuel type.');
end

Q_req_h   = rhs_h + 2502.3*m_evap_h;
Q_total_h = Q_req_h / max(1 - lossFrac, eps);

m_fuel_h = Q_total_h / max(eta_burner*LHV_kJkg, eps);

fuel_per_wetSolid = m_fuel_h / max(Ms_h, eps);
fuel_per_drySolid = m_fuel_h / max(Ss_h, eps);

Q_util_h = rhs_h + 2502.3*m_evap_h;
Q_sup_h  = m_fuel_h * LHV_kJkg * eta_burner;
eta_th   = 100 * Q_util_h / max(Q_sup_h, eps);

fuel = struct();
fuel.fuelType = fuelType;
fuel.LHV_kJkg = LHV_kJkg;
fuel.eta_burner = eta_burner;
fuel.Q_req_h = Q_req_h;
fuel.Q_total_h = Q_total_h;
fuel.m_fuel_h = m_fuel_h;
fuel.fuel_per_wetSolid = fuel_per_wetSolid;
fuel.fuel_per_drySolid = fuel_per_drySolid;
fuel.eta_th = eta_th;

end
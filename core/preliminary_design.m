function design = preliminary_design(Ms, den_s, f_fill, Di, Mgas, rho_g, ...
    Vp, phiD0, Ss_h, Ks, Cps, Cpg, K_gas, mu_g, Yin)

Asec = pi*(Di^2)/4;
Ri   = Di/2;

v_gas = Mgas/(rho_g*Asec);
RPM   = 60*Vp/(pi*Di);

Ss = Ss_h/3600;
rho_s   = den_s;
N_rev_s = RPM/60;

phi = phiD0;
if phi > 1
    phi = phi/100;
end

S_design = 0.3344*Ss /(max(phi,eps)*rho_s*(max(N_rev_s,eps)^0.9)*max(Di,eps));
Ss_design = 100*S_design;
theta_rad = atan(S_design);

Vs = Ms / max(den_s*f_fill*Asec, eps);

betta = ((3*Ms)/(Vs*den_s*(Ri^2)))^(1/3);

Lwcs  = 2*Ri*betta;
Lwncs = Ri*(2*pi - 2*betta);
Lss   = 2*Ri*sin(betta);

De   = (0.5*Di*((2*pi) - 2*betta + sin(2*betta))) / (pi - betta + sin(betta));
Asnc = (pi*(De^2))/4;

nu_g    = mu_g / rho_g;
alpha_g = K_gas / (rho_g*Cpg);
Pr_g    = nu_g / alpha_g;

Re_gas = (v_gas*Di)/nu_g;

if  Re_gas >= 0.4 && Re_gas < 4
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

Nu_cwo = Cc*(Re_gas^mc)*Pr_g^(1/3);
ho     = (Nu_cwo*K_gas)/(Di);

flux_m = Mgas/Asnc;

hgw = Cc*Pr_g^(1/3)*Re_gas^mc;
hgs = hgw*0.10;

Ab  = Ks/(den_s*Cps);
hsw = 0.116*Ks*(((RPM/60)*(Ri^2)*2*betta )/Ab)^0.3 / Di;

design = struct();
design.Asec = Asec;
design.Ri = Ri;
design.v_gas = v_gas;
design.RPM = RPM;
design.phi = phi;
design.S_design = S_design;
design.Ss_design = Ss_design;
design.theta_rad = theta_rad;
design.Vs = Vs;
design.betta = betta;
design.Lwcs = Lwcs;
design.Lwncs = Lwncs;
design.Lss = Lss;
design.De = De;
design.Asnc = Asnc;
design.nu_g = nu_g;
design.alpha_g = alpha_g;
design.Pr_g = Pr_g;
design.Re_gas = Re_gas;
design.Cc = Cc;
design.mc = mc;
design.Nu_cwo = Nu_cwo;
design.ho = ho;
design.flux_m = flux_m;
design.hgw = hgw;
design.hgs = hgs;
design.Ab = Ab;
design.hsw = hsw;

end
function ax = psychrometricChart_VALCON_simplified(varargin)
%PSYCHROMETRICCHART_VALCON_SIMPLIFIED  Carta psicrométrica estilo "VALCON" (simplificada).
%
%   ax = psychrometricChart_VALCON_simplified()
%   ax = psychrometricChart_VALCON_simplified('Name',Value,...)
%
%   + Opcional: representar un punto de estado (Tdb,RH) o (Tdb,w)
%     'PlotState'      true/false (default false)
%     'StateMode'      'T_RH' | 'T_w' (default 'T_RH')
%     'Tdb_C'          (scalar) dry-bulb (°C)
%     'RH_frac'        (scalar) 0..1 (si StateMode='T_RH')
%     'w_state'        (scalar) kg/kg_da (si StateMode='T_w')
%     'StateLabel'     true/false (default true)
%     'StateMarker'    (default 'o')
%     'StateMarkerSize'(default 8)
%     'StateColor'     (default [1 0 0])  % rojo
%
%   (El resto de parámetros son los mismos que antes.)

%% -------------------------
% Parseo de inputs
% -------------------------
p = inputParser;
p.addParameter('Patm_kPa', 101.325, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.addParameter('Tmin', -10, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('Tmax', 55,  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('wmax', 0.040, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.addParameter('RH_list', [0.1 0.3 0.5 0.7 0.9], @(x)isnumeric(x)&&isvector(x));
p.addParameter('Twb_list', 5:5:35, @(x)isnumeric(x)&&isvector(x));
p.addParameter('h_bg_list', 10:10:140, @(x)isnumeric(x)&&isvector(x));
p.addParameter('Figure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('Axes', [], @(x) isempty(x) || ishghandle(x,'axes'));

% ---- Punto de estado (nuevo) ----
p.addParameter('PlotState', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('StateMode', 'T_RH', @(s)ischar(s)||isstring(s));
p.addParameter('Tdb_C', NaN, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('RH_frac', NaN, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('w_state', NaN, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('StateLabel', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('StateMarker', 'o', @(s)ischar(s)||isstring(s));
p.addParameter('StateMarkerSize', 8, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('StateColor', [1 0 0], @(x)isnumeric(x)&&numel(x)==3);

p.parse(varargin{:});
S = p.Results;

Patm_kPa = S.Patm_kPa;
Tmin     = S.Tmin;
Tmax     = S.Tmax;
wmax     = S.wmax;
RH_list  = S.RH_list(:).';
Twb_list = S.Twb_list(:).';
h_bg_list= S.h_bg_list(:).';

if Tmax <= Tmin
    error('Tmax must be greater than Tmin.');
end

%% -------------------------
% Malla de T
% -------------------------
Tvec = linspace(Tmin, Tmax, 900);

%% -------------------------
% Estilos (un color por familia)
% -------------------------
col_sat = [0 0 0];
col_RH  = [0.10 0.10 0.80];
col_Twb = [0.10 0.55 0.10];
col_hbg = [0.55 0.15 0.70];

lw_sat = 2.0; lw_RH = 1.2; lw_Twb = 1.15; lw_hbg = 0.8;

%% -------------------------
% Cálculos
% -------------------------
w_sat = w_saturation(Tvec, Patm_kPa);

w_RH = nan(numel(RH_list), numel(Tvec));
for i = 1:numel(RH_list)
    wtmp = w_from_T_RH(Tvec, RH_list(i), Patm_kPa);
    wtmp(wtmp > w_sat) = NaN;
    wtmp(wtmp < 0)     = NaN;
    w_RH(i,:) = wtmp;
end

w_hbg = nan(numel(h_bg_list), numel(Tvec));
for k = 1:numel(h_bg_list)
    h0 = h_bg_list(k);
    wtmp = (h0 - 1.006*Tvec) ./ (2501 + 1.86*Tvec);
    wtmp(wtmp < 0)     = NaN;
    wtmp(wtmp > w_sat) = NaN;
    w_hbg(k,:) = wtmp;
end

w_Twb = nan(numel(Twb_list), numel(Tvec));
h_Twb = nan(size(Twb_list));
for j = 1:numel(Twb_list)
    Twb = Twb_list(j);
    w_s_wb = w_saturation(Twb, Patm_kPa);
    h_wb   = h_moist_air(Twb, w_s_wb);
    h_Twb(j) = h_wb;

    wtmp = (h_wb - 1.006*Tvec) ./ (2501 + 1.86*Tvec);
    wtmp(wtmp < 0)     = NaN;
    wtmp(wtmp > w_sat) = NaN;
    w_Twb(j,:) = wtmp;
end

%% -------------------------
% Figura/ejes
% -------------------------
if ~isempty(S.Axes)
    ax = S.Axes;
    axes(ax); %#ok<LAXES>
else
    if S.Figure
        figure('Color','w');
    end
    ax = gca;
end

%% -------------------------
% Plot base
% -------------------------
yyaxis(ax,'right'); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

plot(ax, Tvec, w_sat, '-', 'Color', col_sat, 'LineWidth', lw_sat, ...
    'DisplayName','Saturation (RH 100%)');

for k = 1:numel(h_bg_list)
    hv = 'off'; if k==1, hv='on'; end
    plot(ax, Tvec, w_hbg(k,:), ':', 'Color', col_hbg, 'LineWidth', lw_hbg, ...
        'HandleVisibility', hv);
end

for i = 1:numel(RH_list)
    hv = 'off'; if i==1, hv='on'; end
    plot(ax, Tvec, w_RH(i,:), '--', 'Color', col_RH, 'LineWidth', lw_RH, ...
        'HandleVisibility', hv);
end

for j = 1:numel(Twb_list)
    hv = 'off'; if j==1, hv='on'; end
    plot(ax, Tvec, w_Twb(j,:), '-', 'Color', col_Twb, 'LineWidth', lw_Twb, ...
        'HandleVisibility', hv);
end

xlabel(ax, 'Dry-bulb temperature T_{db} (°C)');
ylabel(ax, 'Humidity ratio w (kg_v/kg_{da})');
xlim(ax, [Tmin Tmax]);
ylim(ax, [0 wmax]);
xticks(ax, Tmin:5:Tmax);
title(ax, sprintf('Psychrometric chart (simplified, P = %.3f kPa)', Patm_kPa));

yyaxis(ax,'left');
set(ax,'YTick',[]);
ylabel(ax,'');
yyaxis(ax,'right');

%% -------------------------
% Rotulación
% -------------------------
Tlab_RH = Tmax - 2;
for i = 1:numel(RH_list)
    wlab = w_from_T_RH(Tlab_RH, RH_list(i), Patm_kPa);
    if isfinite(wlab) && wlab <= wmax
        text(ax, Tlab_RH, wlab, sprintf('  %d%%', round(100*RH_list(i))), ...
            'Color', col_RH, 'VerticalAlignment','middle');
    end
end

for j = 1:numel(Twb_list)
    Twb = Twb_list(j);
    Tlab = min(Tmax-5, Twb + 18);
    wlab = interp1(Tvec, w_Twb(j,:), Tlab, 'linear', NaN);
    if isfinite(wlab) && wlab <= wmax
        text(ax, Tlab, wlab, sprintf('  T_{wb}=%d°C', Twb), ...
            'Color', col_Twb, 'Rotation', -35, 'VerticalAlignment','bottom');
    end
end

dx = 0.8; dy = 0.0012;
for j = 1:numel(Twb_list)
    Twb = Twb_list(j);
    if Twb < Tmin || Twb > Tmax, continue; end
    w_on_sat = w_saturation(Twb, Patm_kPa);
    if ~isfinite(w_on_sat), continue; end
    x_txt = Twb + dx;
    y_txt = w_on_sat + dy;
    if y_txt <= wmax
        text(ax, x_txt, y_txt, sprintf('h=%.0f', h_Twb(j)), ...
            'Color', col_sat, 'VerticalAlignment','bottom');
    end
end

legend(ax, {'Saturation (RH 100%)','Enthalpy lines (background)','Relative humidity','Wet-bulb diagonals'}, ...
    'Location','northwest');

%% -------------------------
% Punto de estado (opcional)
% -------------------------
if S.PlotState
    mode = upper(string(S.StateMode));

    if mode == "T_RH"
        if ~isfinite(S.Tdb_C) || ~isfinite(S.RH_frac)
            error('For StateMode="T_RH", provide Tdb_C and RH_frac.');
        end
        T0 = double(S.Tdb_C);
        RH0 = double(S.RH_frac);
        if RH0 < 0 || RH0 > 1
            error('RH_frac must be in [0,1].');
        end
        w0 = w_from_T_RH(T0, RH0, Patm_kPa);

    elseif mode == "T_W"
        if ~isfinite(S.Tdb_C) || ~isfinite(S.w_state)
            error('For StateMode="T_w", provide Tdb_C and w_state.');
        end
        T0 = double(S.Tdb_C);
        w0 = double(S.w_state);

    else
        error('StateMode must be "T_RH" or "T_w".');
    end

    % Validación básica contra saturación y ejes
    w_sat0 = w_saturation(T0, Patm_kPa);
    if ~isfinite(w_sat0)
        error('State point outside saturation function domain.');
    end
    if w0 > w_sat0 + 1e-9
        error('State point is above saturation (w > w_sat at Tdb).');
    end

    yyaxis(ax,'right'); % asegurar eje correcto
    plot(ax, T0, w0, S.StateMarker, ...
        'MarkerSize', S.StateMarkerSize, ...
        'MarkerFaceColor', S.StateColor, ...
        'MarkerEdgeColor', [0 0 0], ...
        'HandleVisibility','off');

    if S.StateLabel
        if mode == "T_RH"
            txt = sprintf('  T=%.1f°C, RH=%.0f%%', T0, 100*RH0);
        else
            txt = sprintf('  T=%.1f°C, w=%.4f', T0, w0);
        end
        text(ax, T0, w0, txt, 'VerticalAlignment','bottom', 'Color', S.StateColor);
    end
end

end

%% =========================
% Funciones auxiliares
% =========================
function Pws_kPa = p_ws_kPa(T_C)
    Pws_kPa = 0.61078 .* exp(17.2694 .* T_C ./ (T_C + 237.29));
end

function w = w_from_T_RH(T_C, RH, P_kPa)
    Pws = p_ws_kPa(T_C);
    Pw  = RH .* Pws;
    Pw  = min(Pw, 0.9999*P_kPa);
    w   = 0.621945 .* Pw ./ (P_kPa - Pw);
end

function w = w_saturation(T_C, P_kPa)
    Pws = p_ws_kPa(T_C);
    Pws = min(Pws, 0.9999*P_kPa);
    w   = 0.621945 .* Pws ./ (P_kPa - Pws);
end

function h = h_moist_air(T_C, w)
    h = 1.006*T_C + w.*(2501 + 1.86*T_C);
end
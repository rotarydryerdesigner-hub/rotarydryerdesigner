function dlg = showPsychrometricChartDialog(parentFig, varargin)

    dlg = uifigure( ...
        'Name','Psychrometric chart', ...
        'WindowStyle','modal', ...
        'Color','white', ...
        'Resize','on', ...
        'Position', local_centerOnParent(parentFig,[980 620]));

    gl = uigridlayout(dlg,[2 1]);
    gl.RowHeight = {'1x', 44};
    gl.ColumnWidth = {'1x'};
    gl.Padding = [10 10 10 10];

    ax = uiaxes(gl);
    ax.Layout.Row = 1;
    ax.Layout.Column = 1;

    btn = uibutton(gl,'Text','Close', ...
        'ButtonPushedFcn', @(~,~) delete(dlg));
    btn.Layout.Row = 2;
    btn.Layout.Column = 1;

    cla(ax,'reset');

    % Dibuja dentro del UIAxes del diálogo (no abre figure externo)
    psychrometricChart_VALCON_simplified('Axes', ax, 'Figure', false, varargin{:});
end

function pos = local_centerOnParent(parentFig, wh)
    pf = parentFig.Position;  % [x y w h]
    w = wh(1); h = wh(2);
    x = pf(1) + (pf(3)-w)/2;
    y = pf(2) + (pf(4)-h)/2;
    x = max(10,x); y = max(10,y);
    pos = [x y w h];
end
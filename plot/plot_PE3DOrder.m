function [fig] = plot_PE3DOrder(Actual, PE3D)
    PE3DOrder  = PE3D.PE3DOrder ;
    bSample    = PE3D.bSample   ;

    color_facecolor = "#FFFFFF";
    color_label     = "#CCCCCC";
    fig_width       = 800;
    fig_height      = 800;
    position = [(1920-fig_width)/2, (1080-fig_height)/2, fig_width, fig_height];

    fig_width_cm  = 10;
    fig_height_cm = 6;
    position_cm   = [5, 5, fig_width_cm, fig_height_cm];  % 显示位置

    figname = 'PE Order';
    fig = figure('Name', figname, 'Position', position, 'Color', color_facecolor);
    % fig = figure('Name', figname, ...
    %              'Units', 'centimeters', ...
    %              'Position', position_cm, ...
    %              'Color', color_facecolor);


    im_PE3DOrder = zeros(size(bSample));

    if strcmpi(Actual.DimFast, 'PE') % Step through Y (#1) first and then Z (#2)
        for ilin = 1:length(PE3DOrder)
            im_PE3DOrder(PE3DOrder(ilin, 1)+PE3D.IdxCenter_PE, PE3DOrder(ilin, 2)+PE3D.IdxCenter_3D) = ilin+1;
        end
        im_PE3DOrder = im_PE3DOrder' ;
    else
        for ilin = 1:length(PE3DOrder)
            im_PE3DOrder(PE3DOrder(ilin, 2)+PE3D.IdxCenter_3D, PE3DOrder(ilin, 1)+PE3D.IdxCenter_PE) = ilin+1;
        end
    end

    cmap = jet(length(PE3DOrder)+1);
    cmap(1,:) = [0,0,0];
    imshow(im_PE3DOrder);      % 显示矩阵

    title(sprintf('R = %d x %d, nRef = %d x %d', Actual.AccelerationPE, Actual.Acceleration3D, ...
        Actual.nRefLinePE, Actual.nRefLine3D), 'FontWeight', 'bold', 'Color', 'r');

    % impixelinfo;
    colormap(cmap);          % 选择颜色映射
    clim([1, length(PE3DOrder)+1]);       % 设置颜色范围
    cb = colorbar;       
    cb.TickLength = 0;

    % 设置 colorbar 样式
    cb.Box = 'on';              
    cb.EdgeColor = 'r';          
    cb.Label.String = 'Order'; 
    cb.Label.Color = 'r';       
    cb.Color = '#CCCCCC';             
    cb.FontSize = 10;           


    ax = gca;
    ax.Position = [0. 0.03 0.85 0.9];      
    ax.LooseInset = [0 0 0 0]; 
    set(fig, 'Color', '#000000');  
    set(fig, 'InvertHardcopy', 'off');  
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0, 0, fig_width_cm, fig_height_cm]);
end
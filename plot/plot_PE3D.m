function [fig] = plot_PE3D(Actual, PE3D)

    PE3DOrder  = PE3D.PE3DOrder ;
    bRef       = PE3D.bRef      ;
    bImaAndRef = PE3D.bImaAndRef;

    fig = figure;
    title(sprintf('Ryz = %d x %d, nRefyz = %d x %d', Actual.AccelerationPE, Actual.Acceleration3D, ...
        Actual.nRefLinePE, Actual.nRefLine3D), 'FontWeight', 'bold', 'Color', 'r');
    xlabel('')
    plot(PE3DOrder(:,1), PE3DOrder(:,2), '.') ;
    hold on
    plot(PE3DOrder(bRef(:),1), PE3DOrder(bRef(:),2), 'o') ;
    plot(PE3DOrder(bImaAndRef(:),1), PE3DOrder(bImaAndRef(:),2), '^', 'MarkerSize', 10) ;

    ax = gca;
    % ax.Position = [0. 0.03 0.80 0.9];      
    ax.LooseInset = [0 0 0 0]; 
    set(fig, 'Color', '#FFFFFF');  
    set(fig, 'InvertHardcopy', 'off');  
end
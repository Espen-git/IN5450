function ploting(x,y,maintitle,type)
    figure()

    % normalize
    y = abs(y)/max(abs(y),[],'all');
    
    % Linear plot
    subplot(2,1,1);
    plot(x,y)
    xlabel('DOA [deg]');
    ylabel('Normalized spectrum');
    title('Linear')
    ylim('padded')
    xlim('tight')
    
    % dB plot
    subplot(2,1,2);
    if ~exist('type','var')
        % Defaults to amplitude
        type = "amplitude";
    end
    
    if type == "power"
        plot(x,10*log10(abs(y)))
    elseif type == "amplitude"
        plot(x,20*log10(abs(y)))
    end
    xlabel('DOA [deg]');
    ylabel('Normalized spectrum [dB]');
    title('dB')
    ylim('padded')
    xlim('tight')
    
    sgtitle(maintitle)
end
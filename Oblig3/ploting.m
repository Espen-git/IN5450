function ploting(image, x, y, maintitle, range)
    figure()
    
    max_value = max((abs(image)),[],'all');
    image_normalized = abs(image) ./ max_value;
    max_value = max((abs(image_normalized)),[],'all');
    
    imshow(db(abs(image_normalized)), [(-range) 0], 'XData', x, 'YData', y)
    colormap turbo
    colorbar
    xlabel('X [m]')
    ylabel('Y [m]')
    title("dB, " + maintitle + ", resolution: " + length(x) + " x " + length(y))
    axis on
end
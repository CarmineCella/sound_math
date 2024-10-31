function testFeatures (filename)
    centroidName = sprintf ('%s.centroid.raw', filename);
    centr = rawread (centroidName, 'double');
    %spread = rawread ('iphone-prototype-deskgoldsmiths.wav.spread.raw', 'double');
    %skewness = rawread ('iphone-prototype-deskgoldsmiths.wav.skewness.raw', 'double');

    energyName = sprintf ('%s.energy.raw', filename);
    energy = rawread (energyName, 'double');

    zcrName = sprintf ('%s.zcr.raw', filename);
    zcr = rawread (zcrName, 'double');

    y = wavread (filename);

    figure
    subplot (4, 1, 1);
    plot (y)
    title ('Sound file')

    subplot (4, 1, 2);
    plot (centr)
    title ('Centroid')

    % subplot (6, 1, 3);
    % plot (spread)
    % title ('Spread')
    % 
    % subplot (6, 1, 4);
    % plot (skewness)
    % title ('Skewness')

    subplot (4, 1, 3);
    plot (zcr)
    title ('ZCR')


    subplot (4, 1, 4);
    plot (energy)
    title ('Energy')


end




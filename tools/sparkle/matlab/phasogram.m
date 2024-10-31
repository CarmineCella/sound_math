function phasogram (file, fftSize, hopSize)
    x = wavread (file);
    empty = zeros (fftSize, 1);
    L = [x; empty];
    pin = 1;
    
    pgram = zeros (fftSize, 1);
    while (pin < length (x)) 
        grain = L(pin:pin+fftSize)'; % .* hanningz(fftSize);
        f = fft (grain, fftSize);
        %mag = abs (f);
        phi = angle (f);
        
        pgram = [ pgram  phi' ];
        pin = pin + hopSize;
    end

    surf (pgram);
end

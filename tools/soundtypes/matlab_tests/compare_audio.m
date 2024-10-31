x = wavread ('cl_low.wav');
or = wavread ('instruments/clarinet.A3.wav');

figure

subplot (2, 1, 1)
[S,F,T,P] = spectrogram (or, 4096, 256, 4096);
surf (T,F,10*log10(P),'edgecolor','none'); axis tight; 
view (0,90);
xlabel ('Time (Seconds)'); ylabel('Hz');
title ('Original sample')

subplot (2, 1, 2)
[S,F,T,P] = spectrogram (x, 4096, 256, 4096);
surf (T,F,10*log10(P),'edgecolor','none'); axis tight; 
view (0,90);
xlabel ('Time (Seconds)'); ylabel('Hz');
title ('Compressed sample')

fname = 'la.wav';
n0 = 10000; %start index
N = 4096; %block length
Nfft = 4096; % FFT length
p = 50; %prediction order
n1 = n0 + N - 1; %end index

[xin, Fs] = wavread (fname, [n0 n1]);
x = xin (:, 1)'; % row vecor of left channel
win = hamming (N)'; % window for input block

[a, g] = lpc (x .* win, p); % calculate LPC coeffs
a
g

A =  g ./ (abs (fft (a, Nfft)));
s =  (abs (fft (x .* win, Nfft)));
H = (A ./ max (A) .* max (s)) + 1;

plot (H(1:end/2) - 1,'r')
hold on
plot (s(1:end/2))

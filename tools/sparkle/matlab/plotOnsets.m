function plotOnsets (wavFile, hop)

[sig, sr] = wavread (wavFile);

markers = rawread ('onsets.raw', 'int');
hfc = rawread ('hfc.raw', 'double');

lay = NaN (length (sig'), 1);
for k = 1 : length (markers)
    lay(k * hop) = markers (k);
end

subplot (2, 1, 1)
plot (sig)
hold on
stem (lay, 'r')
subplot (2, 1, 2)
plot (hfc)
hold on
stem (markers .* max (hfc), 'r')

title ('Signal with time markers');


% eof


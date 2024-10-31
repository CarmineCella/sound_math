function plotSpecEnv (plotTitle, sr)


y = rawread ('orig.raw', 'double');
s = rawread ('transl.raw', 'double');
%p = rawread ('peaksC.raw', 'double');
c = rawread ('ceps.raw', 'double');
cx = rawread ('corr.raw', 'double');

%c = c ./ 10;

ftlen = length (y);
fbin = sr / ftlen;
freqs = [fbin : fbin : sr/2];

figure
plot (freqs, y(1:end/2))
hold on
plot (freqs, s(1:end/2), 'r');
plot (freqs, c(1:end/2), 'k', 'LineWidth', 2);
plot (freqs, cx(1:end/2), 'g');
%plot(freqs,y(1:end/2),'--rs','LineWidth',2,...
 %               'MarkerEdgeColor','k',...
  %              'MarkerFaceColor','g',...
   %             'MarkerSize',10)
grid on

%markers = nan (1, ftlen / 2);
%markers (floor (p) + 1) = y (floor (p) + 1);
%plot (freqs, markers, 'ok', 'LineWidth', 3)

title (plotTitle)
legend ('Magnitude spectrum', 'Translated spectrum', 'Cepstrum');
xlabel ('Hz')
ylabel ('Amp')

axis tight

% eof


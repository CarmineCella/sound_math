% Y = RAWREAD(X, P) Read raw data from a file X with the specified
%                   precision P.
%                   NB: It assumes native byte ordering.

function y = rawread (name, precision)
	prec = '';
	if nargin < 2
		fprintf ('warning: double precision assumed\n');
		prec = 'double';
	else
		prec = precision;
	end
	h = fopen (name, 'r');
	
	if h == -1
		fprintf ('error: cannot open the specified file: %s\n', name);
		y = 0;
		return;
	end
	
	y = fread (h, Inf, prec);
	len = ftell (h);
	fprintf ('%s: %g byte(s) read as %s\n', name, len, prec);
	
	fclose (h);
return;

% EOF


function [a1,a2,a3] = info(sys)

VS=.02; # Criterio 2% valor final
  
[y, t] = step (sys);
[h, w] = size (y);

c=0;

% Obtenemos primer maximo
[ys, i]=max(y);
ts=t(i);
ypico=ys;

##-- Aidan Plenert -- octave-control:stepinfo

% Divide the signal into initial and final segments
% divided on the max point
yhead = y(1:i); % Everything before the peak
ytail = y(i:h); % Everything after peak
ttail = t(i:h);

[hf, wf] = size (ytail);

% Estimate the final convergent value
% Get average of second half of ytail
yfinal = mean (ytail(max (floor (hf/2), 1):hf));

##---------------------------------------------

i=0;
while (abs(ys-yfinal)>VS*yfinal)
  c=c+i+1;
  clear i;
  ## Se busca el maximo en el segmento
  ## del dominio que excluye al ultimo maximo
  ## no satisfactorio
  [ys, i]=max(y(c:h));
  ts=t(i+c);
endwhile

a1=yfinal;a2=ts;a3=ypico;

endfunction

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
j=0;
c=1;
flagMax=0;
flagMin=0;

% Si es sistema sub-amortiguado
if (ypico>1.02*yfinal)

  while (flagMax*flagMin==0)
    c=c+j+i;
    ## Se busca el maximo en el segmento
    ## del dominio que excluye al ultimo maximo
    ## no satisfactorio
    [ys, i]=max(y(c:h));

    % Testeo aceptabilidad del maximo
    if (abs(ys-yfinal)<=VS*yfinal);
      flagMax=1;
    else
      flagMax=0;
    endif

    % Si dicho maximo se encuentra post un minimo aceptable
    if (flagMax*flagMin==1)
      ys=ysm;
      ts=t(j+c);
      break;
    endif

    % Busco minimo inmediato
    [ysm, j]=min(y(c+i:h));

    % Testeo aceptabilidad del minimo
    if (abs(abs(ysm)-yfinal)<=VS*yfinal)
      flagMin=1;
    else
      flagMin=0;
    endif
    
    ts=t(i+c);
    
  endwhile

 else % si NO es sistema sub-amortiguado

   flagMax=1;
   
  while (flagMax*flagMin==0)
  
    c=c+j;

    % busco el primer valor aceptable 
    while (abs(abs(y(c))-yfinal)>VS*yfinal)
      c++;
    endwhile

    % testeamos el minimo inmediato
    [ys, j]=min(y(c:h));

    % verificamos si cumple rango
    if (abs(abs(ys)-yfinal)<=VS*yfinal)
      flagMin=1;
    else
      flagMin=0;
    endif     
    
  endwhile 
  
  ts=t(c);  
 
endif

a1=yfinal;a2=ts;a3=ypico;

endfunction

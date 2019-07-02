##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      |  Respuesta del sistema   | 
##                        |     Aerogenerador        |
##---------------------------------------------------

function [Y,T]=modeloLineal(po,objetivo,atrasoEncendidoPID,J,B,minW,maxW,fu)
  
##########				% CONFIGURACION ENTORNO

  pkg load optim symbolic control;

  warning ('off','OctSymPy:sym:rationalapprox');

  automatic_replot=1;

##########                              % CONFIURACION SOLUCION

  puntosResolucion=5; % Minima cantidad de puntos a obtener para graficar

  errorAdmisible=5; % Porcentaje error admisible de valor alcanzado de RPM
  
##########				% DECLARACIONES

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

  ultimoT=0;
  
  T=0;

  Y=po;

  SET=1;

  t=nan;

  rpmRad=(2*pi)/60;

  radRpm=(1/rpmRad);

  objetivo=objetivo/atrasoEncendidoPID;
  
########## % SCRIPT
  
  while ((1-po/objetivo)*100>=errorAdmisible)
    
    [k,a,yi,yf]=intervaloLineal(minW,maxW,po,fu);
    
    while (SET==1)
      
      [y,t]=laplaceCI(J,B,k,a,po,t);

      SET=0;
      
      if (objetivo>yf || objetivo<yi)  % si el objetivo fuera de intervalo lineal

	maximos=find(y>yf)-1;

	if (isempty(maximos))

	  maximos=find(y<yi)-1;

	endif

      elseif (y(end)>objetivo)   % si objetivo dentro de intervalo

	maximos=find(y>objetivo)-1;

      endif

      if ((isempty(maximos)))  % si no se alcanza el objetivo

	maximo=size(y,1);

      else

	maximo=maximos(1);

      endif   
      
      if (maximo<puntosResolucion)  % si no se cumple resolucion deseada

	t=(0:diff(t)(1)/10:t(end))';

	SET=1;

      endif
      
      
    endwhile
    
    po=y(maximo);
    
    yVector=y(2:maximo,:);

    tActual=t(2:maximo)+ultimoT;

    ultimoT=t(maximo)+ultimoT;
    
    T=[T;tActual];

    Y=[Y;yVector];

    t=nan;SET=1;
    
    figure(2);subplot(2,2,1);plot(T,Y.*(60/(2*pi)),["--" markStyle(1) color(1) ";Simulacion" ";"]);grid on;

    title ("SIMULACION AEROGENERADOR");set(2,"name","MODELO LINEAL");xlabel("t (seg)");ylabel("RPM");
    
  endwhile

  y=Y(2:end);

  t=T(2:end);
  
endfunction


##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      | Respuesta del Sistema    | 
##                        |    a  Lazo Cerrado       |
##---------------------------------------------------

function [y,t]=respuestaLazoCerrado(po,yi,yf,objetivo,J,B,k,a,FR,pid)

##########				% CONFIGURACION ENTORNO

  pkg load optim symbolic control;

  warning ('off','OctSymPy:sym:rationalapprox');

  automatic_replot=1;

##########                              % CONFIURACION SOLUCION

  puntosResolucion=5; % Minima cantidad de puntos a obtener para graficar
  
##########				% DECLARACIONES

  s=tf('s');

  SET=1;

  ZT=0;
  
########## % SCRIPT


  
  while (SET==1)
				% ZERO INPUT


    
    
    thetha=s*(a+s*J*po)/(s*((B-k)+pid*FR)+J*s^2);

    if (ZT==0)   % Si es la primera vez
    
      [ZI,tLocal]=step(thetha);

    else

      [ZI]=step(thetha,tLocal);

    endif
    

				% ZERO STATE


    thetha=(pid*FR)/(J*s+(B-k)+pid*FR);

    [ZE]=step((objetivo)*thetha,tLocal);

    try                    % Posible BUG funcion step, cuando array TLocal proximo a 1000 calcula hasta un elemento menos
          ZT=ZI+ZE;

    catch

      ZT=ZI(1:size(ZE,1))+ZE;

    end_try_catch
    
	  
	  
    SET=0;

    maximos=find(ZT>yf)-1;

    if (isempty(maximos))

      maximos=find(ZT<yi)-1;

    endif

    
    if ((isempty(maximos)))  % si permanece dentro de intervalo

      maximo=size(ZT,1);

    else

      maximo=maximos(1);

    endif   
    
    if (maximo<puntosResolucion)  % si no se cumple resolucion deseada

      tLocal=(0:diff(tLocal)(1)/10:t(end))';

      SET=1;

    endif

  endwhile
  
  t=tLocal(1:maximo);

  y=ZT(1:maximo);
  
endfunction

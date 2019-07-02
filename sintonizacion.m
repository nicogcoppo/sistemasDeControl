##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      | Sintonizacion PID (A&H)  | 
##                        |                          |
##---------------------------------------------------

function [pid,Kc,Ti,Td,alpha]=sintonizacion(po,J,B,k,a,FR)

##########				% CONFIGURACION ENTORNO
  
  pkg load optim symbolic control;

  warning ('off','OctSymPy:sym:rationalapprox');

  automatic_replot=1;

##########				% DECLARACIONES

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

  s=tf('s');

  rpmRad=(2*pi)/60;

  radRpm=(1/rpmRad);
  
########## % ajustes de sintonizacion

  alpha=0.1; % Constante del PID serie
  
  pasoU=20; % Crecimiento de la excitacion en %

  U=5;   % Excitacion inicial

  RELE=50; % Rango de trabajo del RELE (RPM)

  puntosSolucion=4; % Minima cantidad de puntos a plotear de cada rta parcial

  ciclosPrueba=5; % Cantidad de ciclos de prueba para misma excitacion

  criterioConvergencia=5; % Criterio de convergencia general de la solucion en %

  intentosSintonizacion=4;
  
##########				% ADECUACIONES

  rangoControl=RELE;
  
  poInicial=po;
  
  SP=rangoControl*rpmRad; % Valor de activacion rele
  
  ultimoT=0;

  T=0;

  Uglobal=0;

  Y=po;

  A=0;

  flag=1;

  ySET=0;

  instanteCambio=1;

########### % SCRIPT

  
  
  while(flag!=0 && flag<=intentosSintonizacion)

    for i=1:ciclosPrueba

      while (ySET<puntosSolucion)

	


				% ZERO INPUT
	
	
	thetha=(a*s+po*J*s^2)/(J*(s^2)+(B-k)*s);;

	if (ySET==1)    % Primero octave calcula el vector t
	  [ZI]=step(thetha,t);
	else    
	  [ZI,t]=step(thetha);
	endif
	
				% ZERO STATE
	
	
	thetha=(FR)/(J*s+(B-k));
	
	[ZE]=step((-U*mod(i,2))*thetha,t);

	ZT=ZI+ZE;
	
				% Verificacion del calculo numerico
	
	
	test=abs((abs(ZT)-po));

	ySET=find(test>SP,1);

	if (max(test)<=SP)    % Refinamiento del rango de control

	  SP=SP/2;
	  ySET=1;

	elseif (ySET>=puntosSolucion)   % Refinamiento del punto de activacion del rele

	  maximo=ySET-1;

	elseif (ySET<puntosSolucion)

	  t=(0:diff(t)(1)/10:t(end))';
	  ySET=1;
	  
	else
	  
	  ySET=1;
	  maximo=size(ZT,1);
	  
	endif

      endwhile
      
      

  				% PLOTEO
      

           
      figure(4),subplot(2,2,1);set(4,"name","SINTONIZACION");  

      plot(t,1000*ones(size(t)));
      
      plot(t,ZE*radRpm,["--" markStyle(1) color(1) ";ZE;"]);grid on;

      plot(t,ZI*radRpm,["--" markStyle(1) color(2) ";ZI;"]);grid on;
      
      plot(t,ZT*radRpm,["--" markStyle(3) color(3) ";TOTAL;"]);grid on; 

      title ("RESPUESTA PARCIAL");xlabel("t (seg)");ylabel("RPM");
      
				% RESPUESTA GLOBAL

      
      
      po=ZT(maximo);
      
      yVector=ZT(2:maximo,:);

      tActual=t(2:maximo)+ultimoT;

      ultimoT=t(maximo)+ultimoT;
           
      T=[T;tActual];
      
      Y=[Y;yVector];

           
      subplot(2,2,2);
      
      plot(T,Y*radRpm,["--" markStyle(1) color(1) ";sintonizacion;"]);grid on;

      title ("RESPUESTA GLOBAL");xlabel("t (seg)");ylabel("RPM");
      
      subplot(2,2,3);

      Uglobal=[Uglobal;-U*mod(i,2)*ones(size(yVector))];

      Uglobal(end)=0;
      
      plot(T,Uglobal*radRpm,["--" markStyle(1) color(1) ";rele;"]);grid on;

      title ("EXITACION");xlabel("t (seg)");ylabel("RPM");
      
      ySET=0;
      
    endfor

    
    
    U=U*(1+pasoU/100);             % VERIFICA SI EXISTE SINTONIZACION

    maximo=max(Y(instanteCambio:size(Y))*60/(2*pi));

    minimo=min(Y(instanteCambio:size(Y))*60/(2*pi));
    
    if ((maximo-minimo)<=((1+criterioConvergencia/100)*2*rangoControl))   % Verificacion de la convergencia

      flag=0;

    else
      
      instanteCambio=size(Y,1);

      flag++;
      
    endif
    
    
  endwhile



  ku=@(rele,A) (4*rele)/(pi*A);

  kc=@(ku) (ku/1.7);
  
  ti=@(Pu) (Pu/2);

  td=@(Pu) (Pu/8);
  
  Pu=(T(end)-T(instanteCambio))/3;

  PID=@(Ku,Pu,s,a) (Ku/1.7)*((ti(Pu)*s+1)/(ti(Pu)*s))*(td(Pu)/((s*a*td(Pu))+1));

  pid=PID(ku(abs(U),SP),Pu,s,alpha);

  Kc=kc(ku(abs(U),SP));

  Ti=ti(Pu);

  Td=td(Pu);


				% ZERO INPUT
  
  
  
  
  thetha=s*(a+s*J*po)/(s*((B-k)+pid*FR)+J*s^2);

  [ZI,t]=step(thetha,t);

				% ZERO STATE


  thetha=(pid*FR)/(J*s+(B-k)+pid*FR);

  [ZE]=step((poInicial)*thetha,t);

  ZT=ZI+ZE;

  				% PLOTEO


  subplot(2,2,4);

  plot(t,1000*ones(size(t)));

  plot(t,ZE*radRpm,["--" markStyle(1) color(1) ";ZE;"]);grid on;

  plot(t,ZI*radRpm,["--" markStyle(1) color(2) ";ZI;"]);grid on;

  plot(t,ZT*radRpm,["--" markStyle(3) color(3) ";TOTAL;"]);grid on; 

  title ("RESPUESTA PID");xlabel("t (seg)");ylabel("RPM");


    				% RESPUESTA GLOBAL


  maximo=size(ZT,1);

  yVector=ZT(2:maximo,:);

  tActual=t(2:maximo)+ultimoT;

  ultimoT=t(maximo)+ultimoT;

  T=[T;tActual];

  Y=[Y;yVector];

  subplot(2,2,2);

  plot(T,Y*radRpm,["--" markStyle(1) color(1) ";sintonizacion;"]);grid on;

  title ("SINTONIZACION");xlabel("t (seg)");ylabel("RPM");


endfunction

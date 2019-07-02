##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      |     SCRIPT PRINCIPAL     | 
##                        |                          |
##---------------------------------------------------

##########				% CONFIGURACION ENTORNO

pkg load optim symbolic control;

warning ('off','OctSymPy:sym:rationalapprox');

automatic_replot=1;

clear;

##########				% DECLARACIONES

rpmRad=(2*pi)/60;

radRpm=(1/rpmRad);

##########				% PARAMETROS CARACTERISTICOS

J=0.0097200; % Momento de inercia

B=2.40e-3; % Coeficiente de friccion

rho=1.2057; % Densidad

d=0.72; % Diametro

f=1; % Constante de friccion del sistema de frenado

Vinf=5; % Velocidad de la corriente libre

##########                              % CONFIGURACION DE LA SIMULACION

rpmInicial=158; % Revoluciones por minuto minimas para arranque

acceEstacionaria=0.4; % Aceleracion de regimen estacionario


##########                              % PERFORMANCE DEL PID

encendidoPID=95;% Porcentaje de las RPM objetivo que enciende el PID

MS=15; % Maxima sobreelongacion en %

deltaRangoControl=1.1; % Numero a multiplicar al valor del rele para sintonizar

##########                              % FUNCIONES

function [acce,ultimoT,Y,T,po]=ploteo(y,t,Y,T,ultimoT,radRpm)
    
    markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

    color=["k","r","g","b","m","c","k","r","g","b","m","c"];

				% Variables
    
    
    tActual=t+ultimoT;

    ultimoT=t(end)+ultimoT;

    T=[T;tActual];

    Y=[Y;y];   

    A=gradient(Y(2:end))./diff(T);

    acce=A(end);
    
    po=Y(end);

    figure(3);plot(T,Y.*radRpm,["--" markStyle(1) color(1) ";RESPUESTA GLOBAL" ";"]);grid on;

    title ("SIMULADOR");xlabel("t (seg)");ylabel("RPM");set(3,"name","SIMULADOR");

endfunction

##########				% SCRIPT


try     % MENU DE OPCIONES

  eleccion=menu('SELECCIONE OPERACION A REALIZAR','SIMULAR','CALCULAR PARAMETROS PID');

catch

  eleccion=2;

end_try_catch

po=rpmInicial*rpmRad;

[minW,maxW,fu]=interpoladorLsqr(rho,d,Vinf);

###########

if (eleccion==1)    % SIMULADOR
  

  ultimoT=0;

  T=0;

  Y=po;

  i=0;

  MSactual=2*MS;rangoControl=1;

  atrasoEncendidoPID=(encendidoPID/100);

  acce=acceEstacionaria;


  while (po>=0)

    
    
    objetivo=input("INGRESE RPM A MANTENER...")*rpmRad;

    acce=acceEstacionaria;
    
    if (objetivo>atrasoEncendidoPID*po)

      [y,t]=modeloLineal(po,objetivo,atrasoEncendidoPID,J,B,minW,maxW,fu);
      
      [acce,ultimoT,Y,T,po]=ploteo(y,t,Y,T,ultimoT,radRpm);
      
    endif

    while (abs(acce)>=acceEstacionaria)
      
      [k,a,yi,yf]=intervaloLineal(minW,maxW,po,fu);
      
      [pid,Kc,Ti,Td,alpha]=sintonizacion(po,J,B,k,a,f,rangoControl);

      [y,t]=respuestaLazoCerrado(po,yi,yf,objetivo,J,B,k,a,f,pid);
      
      [acce,ultimoT,Y,T,po]=ploteo(y,t,Y,T,ultimoT,radRpm);
      
    endwhile
    
    
    
  endwhile

#########     % CALCULO COEFICIENTES PID

else
  
  [resultados]=resultadosPID(210*rpmRad,700*rpmRad,J,B,f,fu)
  
endif


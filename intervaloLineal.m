##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      |  Generador de intervalos | 
##                        |   lineales segun criterio|
##---------------------------------------------------

function [k,alfa,yI,yF]=intervaloLineal(minW,maxW,po,fu)

  pkg load optim symbolic;

  warning ('off','OctSymPy:sym:rationalapprox');

  
##########				% DECLARACIONES

  syms z x;

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

  i=0;

  errorRel=0;
  
  rpmRad=(2*pi)/60;

  radRpm=(1/rpmRad);

  
##########				% CONFIGURACIONES DEL USUARIO

  errorLinealTolerable=1; % Maximo error de linealizacion tolerable en porcentaje

  resolucionIntervaloLineal=5; % Resolucion de linealizacion en RPM

  pInter=5; % Presicion de la integracion numerica en partes

##########				% SCRIPT

  resolucionIntervaloLinealW=resolucionIntervaloLineal*rpmRad;
  
  funcionSym=fu;

  w=linspace(minW,maxW,100);

  maxIteraciones=(maxW-minW)/(2*resolucionIntervaloLineal);


				% Aproximacion de taylor
  
  derivada=diff(funcionSym,z);

  taylorSym=funcionSym+derivada*(x-z);

  taylor=function_handle(taylorSym);

  funcion=function_handle(funcionSym);

  while (errorRel<=errorLinealTolerable && i<=maxIteraciones)
    
    i++;

    anchoIntervalo=resolucionIntervaloLineal*i;

    limInf=po-anchoIntervalo;

    limSup=po+anchoIntervalo;
    
    if (limInf<minW)
      limInf=minW;
    elseif (limSup>maxW)
      limSup=maxW;
    endif   
      
    espacio=linspace(limInf,limSup,pInter);
    
    areaNoLineal=trapz(espacio,funcion(espacio));

    areaLineal=trapz(espacio,taylor(espacio,po));

    errorRel=abs((areaLineal/areaNoLineal)-1)*100;

  endwhile
  
  anchoIntervalo=resolucionIntervaloLineal*(i-1);

  yI=limInf;

  yF=limSup;
  
  espacio=linspace(yI,yF,10);

  figure(2);
  
  subplot(2,2,4);

  plot(1000*espacio,zeros(length(espacio)));
  
  plot(espacio*radRpm,taylor(espacio,po),["--" markStyle(1) color(1) ";Taylor" ";"]);hold on;grid on;

  plot(espacio*radRpm,funcion(espacio),["--" markStyle(1) color(2) ";Leasqr" ";"]);grid on;

  title ("INTERVALO LINEAL");xlabel("RPM");ylabel("Torque (N m)");
  
  hold off;
  
  coeficientes=coeffs(taylor(x,po),x);

  k=double(coeficientes(1));

  alfa=double(coeficientes(2));
  
endfunction

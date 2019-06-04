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

##########				% DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

##########				% CONFIGURACIONES DEL USUARIO

po=150; % Punto de operacion inicial

pf=875; % Punto de operacion final

J=0.0097200; % Momento de inercia

rho=1.21; % Densidad

d=2.4; % Diametro

Vinf=10; % Velocidad de la corriente libre

##########				% SCRIPT

ultimoT=0;

T=0;

Y=po;

figure(1);subplot(2,2,1);clf;grid on;title ("SIMULACION AEROGENERADOR");  

figure(2);clf;title ("INTERPOLACION CURVA DE TORQUE");

[minRPM,maxRPM,rpmTransicion,fuI,fuF]=interpoladorLsqr(rho,d,Vinf);

while (po<=pf)

  
  [k,a,yi,yf]=intervaloLineal(minRPM,maxRPM,po,rpmTransicion,fuI,fuF);

  [y,t]=laplaceCI(J,k,a,po);
    
  maximo=find(y>yf,1)-1;

  if (y(size(y,1))<yf)  % Si la respuesta se estabiliza y no alcanza el maximo del intervalo
    maximo=size(y,1);
  endif
  

  if (yf==rpmTransicion)
    po=rpmTransicion;
  else
    po=y(maximo);
  endif  
  
  yVector=y(2:maximo,:);

  tActual=t(2:maximo,:)+ultimoT;

  ultimoT=t(maximo)+ultimoT;

  T=[T;tActual];

  Y=[Y;yVector];

  figure(1);subplot(2,2,1);plot(T,Y,["--" markStyle(1) color(1) ";Simulacion" ";"]);grid on;

  title ("SIMULACION AEROGENERADOR");
 

  
endwhile




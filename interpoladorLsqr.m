##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      | Interpolador de Curva    | 
##                        |   de Cp vs. Lambda       |
##---------------------------------------------------

function [minRPM,maxRPM,rpmTransicion,fuI,fuF]=interpoladorLsqr(rho,d,Vinf)

  pkg load symbolic optim;
  
				% DECLARACIONES

  syms z;
  
  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

				% SCRIPT


  data=dlmread('curvaCpTSR.csv',',');

  datum=size(data,1);

  rpm=data(:,1)*(Vinf*60/(2*pi*(d/2)));  

  torque=0.0030447*(data(:,2).*.5*rho*(pi*(d/2)^2)*Vinf^3)./(rpm*2*pi/60);

  minRPM=rpm(1);

  maxRPM=rpm(datum);
  
				% partimos el torque a la mitad

  [maxT,i]=max(torque); % dividimos en el maximo

  rpmTransicion=rpm(i);
  
  torqueI=torque(1:i);rpmI=rpm(1:i);

  torqueF=torque(i:datum);rpmF=rpm(i:datum);

				% Interpolacion Leaser-q/

  leasqrfuncI = @(x, p) p(3) + p(1) * exp (p(2) * x);

  leasqrfuncF = @(x, p) p(3) + p(1)*x + p(2)*sin(x);

  FI = leasqrfuncI;

  FF = leasqrfuncF;

  pinI = [1e-2;1e-2;0];

  pinF = [1e-2;1e-2;maxT]; 

  [fI,pI]=leasqr (rpmI,torqueI, pinI, FI);

  fuI=leasqrfuncI(z,pI);
  
  [fF,pF]=leasqr (rpmF,torqueF, pinF, FF);

  fuF=leasqrfuncF(z,pF);

				% Visualizacion

  figure(2)
  
  hold on;grid on;title ('Regresion no-lineal');
  
  plot(rpmI,torqueI,["--" markStyle(1) color(1) ";Layf Torque" ";"]);

  plot(rpmI,fI,["--" markStyle(1) color(2) ";Leasqr Torque" ";"]);

  plot(rpmF,torqueF,["--" markStyle(2) color(3) ";Layf Torque" ";"]);

  plot(rpmF,fF,["--" markStyle(2) color(4) ";Leasqr Torque" ";"]);grid on;

  hold off;

  ## figure(1);subplot(2,2,2);hold on;grid on;title ('Error');

  ## plot(rpmI,fI./torqueI);

  ## plot(rpmF,fF./torqueF);

  ## hold off

endfunction

  

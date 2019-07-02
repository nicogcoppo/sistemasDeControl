##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      | Interpolador de Curva    | 
##                        |   de Cp vs. Lambda       |
##---------------------------------------------------

function [minW,maxW,fu]=interpoladorLsqr(rho,d,Vinf)
  
  pkg load symbolic optim;
  
				% DECLARACIONES

  syms z;
  
  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

  
  rpmRad=(2*pi)/60;

  radRpm=(1/rpmRad);

				% CONFIGURACION DE LA INTERPOLACION

  nPol=8;
  
				% SCRIPT

  
  
  data=dlmread('curvaCpTSR.csv',',');

  datum=size(data,1);

  w=data(:,1).*(Vinf/(d/2));  

  torque=(data(:,2).*.5*rho*(pi*(d/2)^2)*Vinf^3)./(w);
    
  minW=w(1);

  maxW=w(datum);
  
  % Funcion propuesta
  
  
  
  p=polyfit(w,torque,nPol);

  fu=poly2sym(p,z);
     
				% Visualizacion

  figure(1);clf;set(1,"name","INTERPOLADOR LSQR");subplot(1,2,1);
  
  hold on;grid on;title ('Regresion no-lineal');xlabel("RPM");ylabel("Torque (N m)");
  
  plot(w*radRpm,torque,["--" markStyle(1) color(1) ";Layf Torque" ";"]);

  plot(w*radRpm,polyval(p,w),["--" markStyle(1) color(2) ";Leasqr Torque" ";"]);

  hold off;

  subplot(1,2,2);hold on;grid on;title ('Error');

  xlabel("RPM");ylabel("Error (%)");
  
  plot(w*radRpm,abs(1-polyval(p,w)./torque)*100);

  hold off;

endfunction

  

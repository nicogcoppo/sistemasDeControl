##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      |  Respuesta temporal por  | 
##                        |   el metodo Laplace      |
##---------------------------------------------------


function [y,t]=laplaceCI(J,B,k,a,vo,t)

  
  
  ## CONFIGURACION	
  
  pkg load control

  ## DECLARACIONES

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

  
  rpmRad=(2*pi)/60;

  radRpm=(1/rpmRad);

  
  ## SCRIPT
				% ZERO INPUT

  
  
  s=tf('s');
  
  thetha=(a*s+vo*J*s^2)/(J*(s^2)+(B-k)*s);

  if (isnan(t)==1)       
    [ZI,t]=step(thetha);    
  else
    [ZI]=step(thetha,t);
  endif
  
  y=ZI;

  
  figure(2);
  
  subplot(2,2,2);grid on;
  
  plot(t,ZI*radRpm,["--" markStyle(1) color(1) ";respuesta;"]);grid on;

  title ("RESPUESTA PARCIAL");xlabel("t (seg)");ylabel("RPM");
  

endfunction

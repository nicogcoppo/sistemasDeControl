##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      |  Respuesta temporal por  | 
##                        |   el metodo Laplace      |
##---------------------------------------------------


function [y,t]=laplaceCI(J,k,a,vo)

  ## CONFIGURACION	
  
  pkg load control

  ## DECLARACIONES

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];
  
  ## SCRIPT
				% ZERO INPUT
 
  s=tf('s');
  
  thetha=(a*s+vo*J*s^2)/(J*(s^2)-k*s);
  
  [ZI,t]=step(thetha,18);

  y=ZI;
				% ZERO STATE
  
  ## FT=1/(J*s-k)

  ## [n,d]=numden(FT); ## extraemos numerador y denominador

  ## num=double(coeffs(n)); ## expresamos ambos matricialmente

  ## den=double(coeffs(d)); ## se los convierte en numericos

  ## thetha=tf(num,den); ## planteamos la ft matricialmente para

  ## [ZE]=step(thetha,t);

  figure(1);
  
  subplot(2,2,2);grid on;
  
  plot(t,ZI,["--" markStyle(1) color(1) ";respuesta;"]);grid on;

  title ("RESPUESTA PARCIAL");
  

endfunction

##------* octave *-----------------------------------
##   _____ ____  ______   |    SISTEMAS DE CONTROL   |
##  / ___// __ \/ ____/   |      PROYECTO FINAL      |
##  \__ \/ / / / /        |   Modelo Lineal Aerogen. |
## ___/ / /_/ / /___      |                          |
##/____/_____/\____/      | Compilador de resultados | 
##                        |   de sintonizacion PID   |
##---------------------------------------------------

function [resultados]=resultadosPID(minW,maxW,J,B,FR,fu)

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

  PID=zeros(1,5);    % RPM KC TI TD ALPHA
   
##########                            % SCRIPT

  po=maxW;

  while (po>minW)

    [k,a,yi,yf]=intervaloLineal(minW,maxW,po,fu);
    
    [pid,Kc,Ti,Td,alpha]=sintonizacion(po,J,B,k,a,FR);

    PID=[PID;[po*radRpm Kc Ti Td alpha]];

    po=po-abs(yf-yi)/5;

    figure(5);clf;hold on;

    RPM=PID(2:end,1);
    
    plot(RPM,PID(2:end,2),["--" markStyle(1) color(1) ";KC" ";"]);

    plot(RPM,PID(2:end,3),["--" markStyle(2) color(2) ";TI" ";"]);

    plot(RPM,PID(2:end,4),["--" markStyle(3) color(3) ";TD" ";"]);

    grid on;title ("RESULTADOS COEFICIENTES PID");xlabel("RPM");ylabel("COEFICIENTES");set(5,"name","RESULTADOS PID");

    hold off;
    
  endwhile
  
endfunction

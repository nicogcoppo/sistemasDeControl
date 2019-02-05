######### CONFIGURACION ##########

pkg load control symbolic;

warning ('off','OctSymPy:sym:rationalapprox');

msj=0;save ./octavePidIt msj;

######### DECLARACIONES #######################

FA=sym(zeros(1,1));

FS=sym(zeros(1,1));

PID=sym(zeros(1,1));

FT=sym(zeros(1,1));

solucion=zeros(1,7);

iteracion=zeros(1,7);

syms s c k kp kd ki ;

syms t positive;

######### CARACTERISTICAS MEC. (MKS) ##########

##-- aerogenerador

Ja=2*.2^2; 

Ba=.45;
  
##-- servo-motor

r=0.025;

##-- mordaza de frenado

M=0.5;

Bm=1;

Km=0.5*10/.02;

##-- cable mecanico

Bc=.1;

Kc=.0005*10/.02;  

######### INICIALES/CONTORNO #############

RPM=358*2*pi/60; # velocidad angular inicial

######### FUNCIONES DE TRANSFERENCIA -------
##------- DEL DIAGRAMA DE BLOQUES ########

##-- PID

PID=(864500+302397*s+8822800/s);

##-- AEROGENERADOR

FA=1/(Ja*s+Ba);

##-- SISTEMA DE FRENADO

FS=r*(Kc+s*Bc)/((M*s^2)+s*(Bm+Bc)+Km+Kc);

########## RESPUESTA A UNA PERTURBACION #####

FT=factor(FA/(1+PID*FA*FS)); ## funcion de transf. lazo cerrado

##-- VISUALIZACION DE LA RESPUESTA SIN CONTROL PID

[n,d]=numden(FT); ## extraemos numerador y denominador

num=double(coeffs(n)); ## expresamos ambos matricialmente

den=double(coeffs(d)); ## se los convierte en numericos

thetha=tf(num,den); ## planteamos la ft matricialmente para

step(thetha);

########### ESTIMACION DE PARAMETROS PID #####

##-- OBJETIVO

MP=.05;
TS=.1;
ERROR=.02;

##-- AJUSTES DEL METODO ITERATIVO

% Magnitudes de la var. indep. K 

INICIAL=1;
FINAL=1;
PASO=1;

% Magnitudes de los k dependientes

DIVISIONES=10;
RANGO=10;

##-- PREPARACION DEL PROCESO

k=linspace(INICIAL,FINAL,round((FINAL-INICIAL)/PASO)+1);

kp=linspace(1000000,10000000,DIVISIONES);
kd=linspace(1000000,10000000,DIVISIONES);
ki=linspace(10000000,100000000,DIVISIONES);

##-- ITERACIONES

while size(solucion)(1)<2
  it=0;
  so=0;
  for i=1:size(k,2)
    for j=1:size(kp,2)
      for l=1:size(kd,2)
	for m=1:size(ki,2)
	  
	  try
	    
	    PID=k(i)*(kp(j)+kd(l)*s+ki(m)/s);

			      % Preparamos la funcion de transferencia
	    
	    FT=factor(FA/(1+PID*FA*FS));

	    [n,d]=numden(FT); ## extraemos numerador y denominador

	    num=double(coeffs(n)); ## expresamos ambos matricialmente

	    den=double(coeffs(d)); ## se los convierte en numericos

	    thetha=tf(num,den); ## planteamos la ft matricialmente para

			       % La analisamos y verificamos respuesta
	    
	    [a1,a2,a3]=info(thetha);
	    ys=a1;ts=a2;ypico=a3;
	    
	    
	    if (ypico<=MP && ts<=TS && ys<=ERROR)
	      so++;

              solucion=[solucion;k(i) kp(j) kd(l) ki(m) ypico ts ys];

	      
	      printf "Solucion hallada numero: ",so
	      printf "Sobreelongacion: ",ypico
	      printf "Tiempo de estabilizacion: ",ts
	      printf "Error de estado estacionario: ",ys
	      printf "k= ",k(i)
	      printf "kp= ",kp(j)
	      printf "kd= ",kd(l)
	      printf "ki= ",ki(m)
	      printf "-----------------",so

	      
	      save ./octavePidSol solucion 
	      
	    endif

	    it++;

	    iteracion=[iteracion;k(i) kp(j) kd(l) ki(m) ypico ts ys];

	    mejorIntento=sortrows(iteracion,[6 7])(2,:)([5 6 7]);
	    
	    printf "Iteraciones: ",it
	    printf "Soluciones: ",so
	    printf "Sobreelongacion: ",ypico
	    printf "Tiempo de estabilizacion: ",ts
	    printf "Error de estado estacionario: ",ys
	    printf "k= ",k(i)
	    printf "kp= ",kp(j)
	    printf "kd= ",kd(l)
	    printf "ki= ",ki(m)
	    printf ("Mejor intento -> M=%d ts=%d Err=%d \n",mejorIntento(1),mejorIntento(2),mejorIntento(3))
	    printf "-----------------",it

	    [fi,c]=size(iteracion);
	    msj=iteracion(fi,:);
	    save ./octavePidIt -append msj 

	end_try_catch
	
      endfor
    endfor
  endfor
endfor

% NUEVO DOMINIO DE ANALISIS EN FUNCION DE LA MEJOR OPCION

dominioInicial=sortrows(iteracion,[6 7])(2,:)([2 3 4])*(1-RANGO/100);

dominioFinal=sortrows(iteracion,[6 7])(2,:)([2 3 4])*(1+RANGO/100);

kp=linspace(dominioInicial(1),dominioFinal(1),DIVISIONES);

kd=linspace(dominioInicial(2),dominioFinal(2),DIVISIONES);

ki=linspace(dominioInicial(3),dominioFinal(3),DIVISIONES);

endwhile


  

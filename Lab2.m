clc;
clear all;
close all;

%Datos a utilizar
f_0 = 10; %Frecuencia Señal
fm = 1000; %Frecuencia de muestreo
t=(-f_0:1/fm:f_0);
T=1/f_0;
Ts=1/fm;
%Muestras
muestras=T/Ts;


%Señal Sen
x = 2*f_0*pi*t;
seno = sin(x);
sac = seno./x;
idx = isnan(sac);
sac(idx) = 1;
%Creacion de coseno
r=0;
f_d = f_0*r;
x = 2*f_d*pi*t;
cosen = cos(x);
Den = (1-(4*f_d*t).^2);
cost = cosen./Den;
h = 2*f_0*sac.*cost;%Coseno alzado

r = 0.25;
f_d = f_0*r;
x = 2*f_d*pi*t;
cosen = cos(x);
Den = (1-(4*f_d*t).^2);
cost = cosen./Den;
h1 = 2*f_0*sac.*cost;%Coseno alzado

r = 0.75;
f_d = f_0*r;
x = 2*f_d*pi*t;
cosen = cos(x);
Den = (1-(4*f_d*t).^2);
cost = cosen./Den;
h2 = 2*f_0*sac.*cost;%Coseno alzado

r = 1;
f_d = f_0*r;
x = 2*f_d*pi*t;
cosen = cos(x);
Den = (1-(4*f_d*t).^2);
cost = cosen./Den;
h3 = 2*f_0*sac.*cost;%Coseno alzado

figure(1)
plot(t,h,'b','LineWidth',2);
hold on
plot(t,h1,'m','LineWidth',2);
plot(t,h2,'r','LineWidth',2);
plot(t,h3,'c','LineWidth',2);
legend('roll-off=0','roll-off=0.25','roll-off=0.75','roll-off=1');
grid on
xlim([-7*T 7*T]);%Graficar solo 7 periodos
xlabel('Tiempo')
ylabel('Amplitud')
title('Coseno alzado')


%Ahora sacamos la respuesta en frecuencia
L = length(h);%Largo de datos de tiempo
n = 2^nextpow2(L);
f = f_0*(0:(n/2))/n;

%|H|
fft_h = fft(h,n);
P_h = abs(fft_h/n).^2;

%|H1|
fft_h = fft(h1,n);
P_h1 = abs(fft_h/n).^2;

%|H2|
fft_h = fft(h2,n);
P_h2 = abs(fft_h/n).^2;

%|H3|
fft_h = fft(h3,n);
P_h3 = abs(fft_h/n).^2;



figure(2)
plot(f,P_h(1:n/2+1),'b','LineWidth',2);
hold on
plot(f,P_h1(1:n/2+1),'m','LineWidth',2);
plot(f,P_h2(1:n/2+1),'r','LineWidth',2);
plot(f,P_h3(1:n/2+1),'c','LineWidth',2);
legend('roll-off=0','roll-off=0.25','roll-off=0.75','roll-off=1');
grid on
axis([0 0.25 0 0.0015])
xlabel('Frecuencia')
ylabel('|H|')
title('Respuesta en Frecuencia Coseno alzado')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creacion de coseno alzado para esta parte con r=0.22
fx=f_0/2;%6db de f_0
x = 2*fx*pi*t;
seno = sin(x);
sac = seno./x;
idx = isnan(sac);
sac(idx) = 1;
%Cos
r = 0.22;
f_d = fx*r;
x = 2*f_d*pi*t;
cosen = cos(x);
Den = (1-(4*f_d*t).^2);
cost = cosen./Den;
h4 = 2*fx*sac.*cost;%Coseno alzado r=0.22

sym=10^5;%simbolos
%Creacion de numeros aleatorios
aleatbinaria = randi([0,1], [1,sym]);
%Pasar los valores de 0 a -1
for i=1:length(aleatbinaria)
    if aleatbinaria(i) == 0
        aleatbinaria(i) = -1;
    end
end

t1=[-(length(aleatbinaria)/2)*T:Ts:(length(aleatbinaria)/2)*T];

%Impulsos
Impulsos=upsample(aleatbinaria,muestras);

%Grafico de coseno alzado para BPSK y impulsos
figure(3)
plot(t,h4,'-o');
xlim([-7*T 7*T]);
xlabel('Tiempo'); 
ylabel('Amplitud');

filtro=conv(Impulsos,h4);

t3=[-((length(filtro)/muestras)/2*T):Ts:(length(filtro)/muestras)/2*T];

figure(4);
subplot(1,1,1); 
plot(t1, [Impulsos 0],'-o');
xlim([-20*T 20*T]); 
title('Impulsos'); 
xlabel('Tiempo'); 
ylabel('Amplitud');
figure(5);
subplot(1,1,1); 
plot(t3, [filtro 0],'-o');
title('Coseno alzado filtrado(solo 20 periodos)');
xlim([-20*T 20*T]); 
xlabel('Tiempo');
ylabel('Amplitud');

filtro_reshape = reshape(filtro,2*muestras,[]);
AWNG=awgn(filtro_reshape,15); 

figure(6);
subplot(1,1,1)
plot(filtro_reshape,'b');   
title('Diagrama de Ojo roll-off=0.22');
xlabel('Tiempo')
ylabel('Amplitud')
grid on
figure(7);
subplot(1,1,1);
plot(AWNG,'b');
title('Diagrama de ojo AWNG roll-off=0.22');
xlabel('Tiempo')
ylabel('Amplitud')
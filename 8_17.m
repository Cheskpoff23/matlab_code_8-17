close all, clear, clc


%limites 
xm = -700; xM = 1500 ;
ym = -500; yM = 1000;


p21 = 667.55077098 + 156.00717365i ;
p31 = 1093.88859278 + 888.84591052i ;
O2 = 349.70640626 + 573.71400323i;
O4 = 49.70640626 +573.71400323i; 
a2 = -4.1;
a3 = 0;


angulos = [];
for i = 0:1 
if i == 0 
    
else 
    O2 = O4;
end 

p1 = 0;
R1 = p1 - O2 ;
R2 = R1 + p21;
R3 = R1 + p31;

w1 = atand(imag(R1)/real(R1));
if real(R1) < 0  
    w1 = 180 + w1; 
end 
w2 = atand(imag(R2) / real(R2));
if real(R2) < 0  
    w2 = 180 + w2; 
end 
w3 = atand(imag(R3) / real(R3));
if real(R3) < 0  
    w3 = 180 + w3; 
end 
C1 = norm(R3)*cosd(a2 + w3) - norm(R2)*cosd(a3 + w2);
C2 = norm(R3)*sind(a2 + w3) - norm(R2)*sind(a3 + w2);
C3 = norm(R1)*cosd(a3 + w1) - norm(R3)*cosd(w3);
C4 = -norm(R1)*sind(a3 + w1) + norm(R3)*sind(w3);
C5 = norm(R1)*cosd(a2 + w1) - norm(R2)*cosd(w2);
C6 = -norm(R1)*sind(a2 + w1) + norm(R2)*sind(w2);

A1 = -(C3^2)-(C4^2);
A2 = C3*C6 - C4*C5;
A3 = -C4*C6 - C3*C5;
A4 = C2*C3 + C1*C4;
A5 = C4*C5 - C3*C6;
A6 = C1*C3 - C2*C4;

K1 = A2*A4 + A3*A6;
K2 = A3*A4 + A5*A6;
K3 = ((A1^2)-(A2^2)-(A3^2)-(A4^2)-(A6^2)) / 2;


if i == 0

op_1_b3 = 2*atand((K2 - sqrt((K1^2)+(K2^2)-(K3^2))) / (K1 + K3));
op_2_b3 = 2*atand((K2 + sqrt((K1^2)+(K2^2)-(K3^2))) / (K1 + K3));

if abs(op_2_b3 - a3) < 0.01
    b3 = op_1_b3;
else 
    b3 = op_2_b3;
end 

%----------------gamma---------------%
else 
op_1_b3 = 2*atand((K2 - sqrt((K1^2)+(K2^2)-(K3^2))) / (K1 + K3));
op_2_b3 = 2*atand((K2 + sqrt((K1^2)+(K2^2)-(K3^2))) / (K1 + K3));

if abs(op_2_b3 - a3) < 0.01
    b3 = op_1_b3;
else 
    b3 = op_2_b3;
end 
end 

b2 = atand(-(A3*sind(b3) + A2*cosd(b3) + A4) /  -(A5*sind(b3) + A3*cosd(b3) + A6));

if -(A5*sind(b3) + A3*cosd(b3) + A6) < 0  
    b2 = 180 + b2; 
end
angulos = [angulos ; b2 b3]; %en la segunda parte del if arroja los valores de g2 y g3
end 

b2 = deg2rad(angulos(1)); 
b3 = deg2rad(angulos(3));
g2 = deg2rad(angulos (2)); 
g3 = deg2rad(angulos(4));
%----pequeña correccion con los angulos alfa 2 y 3 de entrada
a2 = deg2rad(a2);
a3 = deg2rad(a3);


%% construcción de la matriz y solución [W Z U S]

A = [(exp(1i*b2)-1) (exp(1i*a2)-1) ; (exp(1i*b3)-1) (exp(1i*a3)-1)];
A2 = [(exp(1i*g2)-1) (exp(1i*a2)-1) ; (exp(1i*g3)-1) (exp(1i*a3)-1)];

B = [p21; p31];
%B = [p21*exp(1i*d2); p31*exp(1i*d3)];
sol = cat(1, A\B, A2\B);

W = sol(1); Z = sol(2);
U = sol(3); S = sol(4);

%%  [W Z U S] EXTRAIDOS DIRECATMENTE DE LAS COTAS

%W = -718 + 291i;
%V = 340 - 60i;
%U = -922 + 231i;
%S = 740 + 60i;
%Z = 1080 + 0i;


% ubicación puntos de los pivotes
O2 = 0 - Z - W;
O4 = -S - U;
% construcción de las barras restantes
V = Z - S;
G = W + V - U;


%% animación 
%magnitud de las barras

a = norm(W);
b = norm(V);
c = norm(U);
d = norm(G);
e = norm(Z);
f = norm(S);

%% gráfica estática
%figure() %simepre fuera del for :)
% barras

grid on
gM4b = [O2 O2+W O2+W+V O2+W+V-U];
plot(gM4b,'-ro','LineWidth',3)
hold on
gZ = [O2+W O2+W+Z];
plot(gZ,'-g','LineWidth',3)
gS = [O2+W+V O2+W+V+S];
plot(gS,'-g','LineWidth',3)
gG = [O2 O2+G];
plot(gG,'-ko','LineWidth',3)

xlim([xm xM])
ylim([ym yM])

% puntos de precisión
%plot(delta2,'bo','LineWidth',1)
%plot(delta3,'bo','LineWidth',1)
plot(0,0,'bo','LineWidth',1)
grid on
hold off

pause (5) %para visualizar la gráfica estática

vV = [real(V) imag(V)];
vZ = [real(Z) imag(Z)];

Pp = (vV(1) * vZ(1)) + (vV(2) * vZ(2));
R3_5 = acosd( Pp /(norm(vV) * norm(vZ)));

%angulos  thetha 1 y thetha 2

vector_tierra = O4 - O2;

R1 = atand(imag(vector_tierra) / real(vector_tierra));  % de la tierra
if real(vector_tierra) < 0 
   R1 = 180 + R1; 
end 
%-----------------T2
R2_inicial = atand(imag(W) / real (W));% de la barra W
if real(W) < 0 
    R2_inicial = R2_inicial + 180;
end 
%---------------T3
R3_aprox = atand(imag(V) / real (V));% de la barra V ;
if real(V) < 0  
    R3_aprox = 180 + R3_aprox;
end 
%--------T4
R4_aprox = atand(imag(U) / real (U));% de la barra U ;
if real(U) < 0 
    R4_aprox = 180 + R4_aprox;
end

matriz_analisis = [];
matriz_graf = [];
matriz_graf_acopladora = [];
rango = [];
transmision = [];
theta2 = [];
for R2 = R2_inicial:1:(R2_inicial+93) %R2_inicial-360:1:R2_inicial R2_inicial:-1:-R2_inicial %-69.96509800916945
  
    cla
    %% Newton Raphson [W V U G]

    
    
    R3 = R3_aprox;
    R4 = R4_aprox;
       x0 = [R3 ; R4];

    

    % Definir la tolerancia
    tolerancia = 1e-9;
    
    % Número máximo de iteraciones
    max_iter = 10000;
    
    % Inicializar el contador de iteraciones
    iter = 0;
    
    while iter < max_iter
        f1 = a*cosd(R2) + b*cosd(x0(1)) - c*cosd(x0(2)) - d*cosd(R1);
        f2 = a*sind(R2) + b*sind(x0(1)) - c*sind(x0(2)) - d*sind(R1);
    
        df1_dR3 = -b*sind(x0(1));
        df1_dR4 = c*sind(x0(2));
        df2_dR3 = b*cosd(x0(1));
        df2_dR4 = -c*cosd(x0(2));
    
        % Calcular los valores de las funciones y el jacobiano en el punto actual
        F = [f1 ; f2];
        J = [df1_dR3 df1_dR4 ; df2_dR3 df2_dR4];
        
        % Calcular la nueva aproximación usando el método de Newton-Raphson
        x1 = x0 - (J\F);
        
        % Verificar la convergencia
        if norm(J\F) < tolerancia
            xx = 1;
            %disp("Convergió");
            rango = [rango ; R2];
            break;
         end
        
        % Actualizar el punto de aproximación
        x0 = x1;
        
        % Incrementar el contador de iteraciones
        iter = iter + 1;
    end
    R3 = x1(1);
    R4 = x1(2);
while R3 < 0
    R3 = R3 + 360;
end
while R3 > 360 
    R3 = R3 - 360;
end 
while R4 < 0
    R4 = R4 + 360;
end
while R4 > 360 
    R4 = R4 - 360;
end 
         
    if iter == max_iter
        xx = 0;

      %disp('El método de Newton-Raphson no convergió.');
      
    end
         R3_aprox = R3;
         R4_aprox = R4;
    %% el resto :) [V S Z]
%% menos R3_5 si el gorro esta abajo, mas si está arriba.
   R5 = R3 - R3_5;   

    while R5 < 0
    R5 = R5 + 360;
    end
    while R5 > 360 
    R5 = R5 - 360;
    end 


   
    %% animacion 
    if xx == 1 
    hold on

    ZG=d*(cosd(R1)+1i*sind(R1)); 
    ZW=a*(cosd(R2)+1i*sind(R2));
    ZV=b*(cosd(R3)+1i*sind(R3));
    ZU=c*(cosd(R4)+1i*sind(R4));
    ZZ=e*(cosd(R5)+1i*sind(R5)); %ZZ=e*(cosd(R5)+1i*sind(R5)); %ZZ=e*(cosd(x1_2(1))+1i*sind(x1_2(1)));
    ZS = ZZ - ZV;%ZS=f*(cosd(R6)+1i*sind(R6));
  %-----para hallar el angulo T6 en cada posicion  -----%

    R6 = atand(imag(ZS) / real(ZS));
    if real(ZS) < 0 
        R6 = R6 + 180;
    end

    matriz_analisis = [matriz_analisis; R2 R3  R4  R5 R6 ];

    
   
    %----------grafica del punto P---------%
    p = O2+ZW+ZZ;
    matriz_graf = [matriz_graf; p ];
    plot (matriz_graf);
    %-----------grafica puntomedio de la barra acopladora------%
    %p_acopladora = O2+ZW+(ZV/2);
    %matriz_graf_acopladora = [matriz_graf_acopladora; p_acopladora ];
    %plot (matriz_graf_acopladora);

    %---para que me diga cual es la posicion ---3%

    if abs(p - p31) < 0.01
        AAAposicion_3 = R2;

        while AAAposicion_3 < 0
   AAAposicion_3 = AAAposicion_3 + 360;
    end
    while AAAposicion_3 > 360 
   AAAposicion_3 = AAAposicion_3 - 360;
    end 

    end 

    tierra = [O2 O4];
    plot(tierra,'-ko','LineWidth',3)
    
    mecanismo = [O2 O2+ZW O2+ZW+ZV O2+ZW+ZV-ZU];
    plot(mecanismo,'-ro','LineWidth',3)

    acoplador = [O2+ZW O2+ZW+ZZ O2+ZW+ZZ-ZS];
    plot(acoplador,'-go','LineWidth',3)

    %PUNTOS DE PRESICION 
    plot(p21,'bo','LineWidth',1)
    plot(p31,'bo','LineWidth',1)
    plot(0,0,'bo','LineWidth',1)

    
    xlim([xm xM])
    ylim([ym yM])


    grid on
    pause(0.2)

    %ángulo de transmición y theta 2
    pp = real(ZV)*real(ZU) + imag(ZV)*imag(ZU) ;
    mag = norm(ZV)*norm(ZU);
    Rt = acosd(pp / mag);

    transmision = [transmision ; Rt ];
    theta2 = [theta2 ; R2];
    end 
   
end
figure()
plot(theta2,transmision)
title('Tetha 2 vs ángulo de transmisión')
rango_angulo_transmision = [min(transmision) , max(transmision)]

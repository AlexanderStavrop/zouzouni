clc
clear;
close all;
format long;
tic

%% ~~~~~~~~~~~~~~~~~~~~~~~ particle much denser than the surrounding fluid ~~~~~~~~~~~~~~~~~~~~~~~

St=0.001; % [0.001 0.5 10]
rc=1; %kg/m3 ,density of carrier phase (flow)
rd=10; % kg/m3 density of dispersed phase
Dp=10^(-3);
rp=Dp/2;
mp=rd*pi/6*Dp^3;
mf=1;
Uo=1;
mc=10^(-3); % dynamic viscosity of carrier phase (air)
g=9.81; %gravity
a=0.63661977236; %*(pi/2);
R=mf/(mp+mf/2);

%% Topology

N=20;                   % Number of particles in each row/line
time = 1;

XMIN=-1;
XMAX=3;
YMIN=0;
YMAX=4;

dx=(XMAX-XMIN)/N;
dy=(YMAX-YMIN)/N;


% Initializing needed matrices
x1 = zeros(N,N);
x2 = zeros(N,N);

u1 = zeros(N,N);
u2 = zeros(N,N);

V1 = zeros(N,N);
V2 = zeros(N,N);


%% ~~~~~~~~~~~~~~~~~~~~~~ Equations of Motion ~~~~~~~~~~~~~~~~~~~~~~
for i=1:N
    for j=1:N
        % Calculating x1
        x1(i,j)= XMIN+j*dx;
        
        % Calculating X2
        x2(i,j)= YMIN+i*dy;

        % Calculating u1 and u2
        u1(i,j) = Uo*cos(x1(i,j)/a)*cos(x2(i,j)/a); % fluid velocity in x direction
        u2(i,j) = Uo*sin(x1(i,j)/a)*sin(x2(i,j)/a); % fluid velocity in y direction
    
        % Calculating V1 and V2
        V1(i,j) = u1(i,j);% particle velocity in x direction
        V2(i,j) = u2(i,j); % particle velocity in y direction
    end
end


%%

dt=0.001;
Rer1 = zeros(N,N);
Rer2 = zeros(N,N);

f1 = zeros(N,N);
f2 = zeros(N,N);

steady_state_x = zeros(N,N);
external_forces_x = zeros(N,N);
virual_mass_x = zeros(N,N);

steady_state_y = zeros(N,N);
external_forces_y = zeros(N,N);
virual_mass_y = zeros(N,N);

dV1 = zeros(N,N);
dV2 = zeros(N,N);

V1_new = zeros(N,N);
V2_new = zeros(N,N);

x1_new = zeros(N,N);
x2_new = zeros(N,N);


for t=1:time

    for i=1:N
        for j=1:N
            value1 = abs(u1(i,j)-V1(i,j))*Dp*rc/mc;
            value2 = abs(u2(i,j)-V2(i,j))*Dp*rc/mc;
            
            % First matrix
            if value1==0
                Rer1(i,j)=10^(-7);
            else 
                Rer1(i,j)=value1;
            end
            f1(i,j)= 1;%+0.15*Rer1(i,j).^0.687;

            % Second matrix
            if value2==0
                Rer2(i,j)=10^(-7);
            else 
                Rer1(i,j)=value1;
            end
            f2(i,j)= 1;

            % Calculating dv1
            tv = rd*Dp^2/(18*mc);
            steady_state_x(i,j) = 1/tv*(Uo*cos(x1(i,j)/a)*cos(x2(i,j)/a)-V1(i,j));
            external_forces_x(i,j) = -rc/rd*(Uo^2/a*cos(x1(i,j)/a)*sin(x1(i,j)/a));   
            virual_mass_x(i,j) = rc/rd*(-V1(i,j)/2*Uo/a*sin(x1(i,j)/a)*cos(x2(i,j)/a)-V2(i,j)/2*Uo/a*cos(x1(i,j)/a)*sin(x2(i,j)/a));

            dV1(i,j) = (1/(1+rc/rd))*(steady_state_x(i,j) + external_forces_x(i,j) + virual_mass_x(i,j))*dt;


            % Calculating dv2
            gravity=-g;
            steady_state_y(i,j) = 1/tv*(Uo*sin(x1(i,j)/a)*sin(x2(i,j)/a)-V2(i,j)/2);
            external_forces_y(i,j) = rc/rd*(Uo^2/a*cos(x2(i,j)/a)*sin(x2(i,j)/a)+g);
            virual_mass_y(i,j) = rc/rd*(V1(i,j)/2*Uo/a*cos(x1(i,j)/a)*sin(x2(i,j)/a) + V2(i,j)/2*Uo/a*sin(x1(i,j)/a)*cos(x2(i,j)/a));    

            dV2(i,j) = (1/(1+rc/rd))*(gravity+steady_state_y(i,j) + external_forces_y(i,j) + virual_mass_y(i,j))*dt;
        
            
            V1_new(i,j) = V1(i,j) + dV1(i,j);
            V2_new(i,j) = V2(i,j) + dV2(i,j);

            x1_new(i,j) = x1(i,j) + (V1(i,j) + V1_new(i,j))*dt/2;
            x2_new(i,j) = x2(i,j) + (V2(i,j) + V2_new(i,j))*dt/2;
    
        end
    end
    

    for i=1:N
        for j=1:N
            if x1_new(i,j)>XMAX
                x1_new(i,j)=XMIN + x1_new(i,j)-XMAX;
            elseif x1_new(i,j)<XMIN
                x1_new(i,j)=XMAX-abs(x1_new(i,j)-XMIN);
            end

            if x2_new(i,j)>YMAX
                x2_new(i,j)=YMIN + x2_new(i,j)-YMAX;
            elseif x2_new(i,j)<YMIN 
                x2_new(i,j)=YMAX-abs(x2_new(i,j)-YMIN);              
            end
        end
    end
  
    x1=x1_new;
    x2=x2_new;
    V1=V1_new;
    V2=V2_new;
    u1=Uo*cos(x1./a).*cos(x2./a); % fluid velocity in x direction
    u2=Uo*sin(x1./a).*sin(x2./a); % fluid velocity in y direction
  

    figure();
   
    for i=1:N
        for j=1:N
            hold on
            scatter(x1(i,j),x2(i,j),5,'k','o'); %Topology
            axis([XMIN XMAX YMIN YMAX]);
        end
    end
%     frame(t)=getframe(gcf);
%     disp("done with " + t)
end

toc
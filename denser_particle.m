clc
clear all
format long
                        %~~~~~~~~ particle much denser than the surrounding fluid ~~~~~~~~~~~~~~~~~~~~~~~        

R=1*10^(-3); % R=(mf/(mp+1/2mf)
%Q=0.01;
%B=0.5;
%A=2; %A=1/Stk
St=0.001% [0.001 0.5 10] 

rc=1; %kg/m3 ,density of carrier phase (flow)
rd=rc*(1-R/2)/R; % kg/m3 density of dispersed phase
Dp=10^(-6);
rp=Dp/2;
mp=rd*pi/3*Dp^3;
mf=rc*pi/3*Dp^3;
mc=1.78*10^(-5); % dynamic viscosity of carrier phase (air)
g=9.81; %gravity
%Uo=(mp-mf)*g/(6*pi*Dp/2*mc*Q);
L=0.05;
Uo=18*mc*St*L/(rd*Dp^2);
a=L*pi;
%B=(6*pi*Dp/2*mc*L)/((mp+mf/2)*Uo);

N=10; %number of particles in each row/line
XMIN=-L;
XMAX=L;
YMIN=-L;
YMAX=L;

dx=(XMAX-XMIN)/N;
dy=(YMAX-YMIN)/N;

%Topology

for i=1:N
    for j=1:N
        x1(i,j)= XMIN+j*dx;
    end
end
for i=1:N
    for j=1:N
        x2(i,j)= YMIN+i*dy;
    end
end
%~~~~~~~~~~~~~~~~ Graphic representation of Topology ~~~~~~~~~~~~~~~~
%x2=x1';
%hold on
% for i=1:N
%     for j=1:N
%            
%         hold on;
%         scatter(x1(i,j),x2(i,j),5,'k','o'); %Topology
%     end
%     axis([XMIN-0.5 XMAX+0.5 YMIN-0.5 YMAX+0.5]);
% end
% hold off


% XMAX=XMAX/L;
% XMIN=XMIN/L;
% YMAX=YMAX/L;
% YMIN=YMIN/L;
% x1=x1/L;
% x2=x2/L;
            %~~~~~~~~~~~~~~~~~~~~~ Equations of Motion~~~~~~~~~~~~~~~~~~~~~~

u1=Uo*cos(x1./a).*cos(x2./a); % fluid velocity in x direction
u2=Uo*sin(x1./a).*sin(x2./a); % fluid velocity in y direction
V1=Uo*cos(x1./a).*cos(x2./a); % particle velocity in x direction
V2=Uo*sin(x1./a).*sin(x2./a); % particle velocity in y direction

dt=0.001;
Rer1=zeros(N,N);
Rer2=zeros(N,N);

% v=VideoWriter('plot2.avi');
% open(v);

for t=1:2
    for i=1:N
        for j=1:N
            Rer1(i,j)=abs(u1(i,j)-V1(i,j))*Dp*rc/mc;
            Rer2(i,j)=abs(u2(i,j)-V2(i,j))*Dp*rc/mc;
        end
    end

    for i=1:N
        if Rer1(i,j)==0
           Rer1(i,j)=10^(-7);
        end
        if Rer2(i,j)==0
           Rer2(i,j)=10^(-7);
        end
    end
        f1=1+0.15*Rer1.^0.687;
        f2=1+0.15*Rer2.^0.687;
    
        if (Rer1>800) 
            22222222 % indicative of an error
        end
        if (Rer2>800) 
            22220000 %indicative of an error
        end

    dV1=(-6*pi*rp*mc*f1.*V1./(mp+mf/2)+(6*pi*rp*mc*Uo*f1./(mp+mf/2))*cos(x1./a).*cos(x2./a)-R*(Uo^2/a)*cos(x1./a).*sin(x1./a)-1/2*R*(Uo/a)*(V1.*sin(x1./a).*cos(x2./a)+V2.*cos(x1./a).*sin(x2./a))).*dt;
    dV2=(((mf-mp)/(mp+1/2*mf))*g-6*pi*rp*mc*f2.*V2./(mp+1/2*mf)+6*pi*rp*mc*f2.*Uo*sin(x1./a).*sin(x2./a)/(mp+mf/2)+R*2*Uo^2/a*cos(x2./a).*sin(x2./a)+Uo/(2*a)*R*(V1.*cos(x1./a).*sin(x2./a)+V2.*sin(x1./a).*cos(x2./a))).*dt;
    
    V1_new=V1+dV1;
    V2_new=V2+dV2;

    x1_new=x1+(V1+V1_new).*dt/2;
    x2_new=x2+(V2+V2_new).*dt/2;


    for i=1:N
        for j=1:N
            if x1_new(i,j)>XMAX
                x1_new(i,j)=XMIN + x1_new(i,j)-XMAX;
                
               1;
            end
            if x1_new(i,j)<XMIN
                x1_new(i,j)=XMAX-abs(x1_new(i,j)-XMIN);
                
                2;
            end
            if x2_new(i,j)>YMAX
                x2_new(i,j)=YMIN + x2_new(i,j)-YMAX;
                3;
            end
            if x2_new(i,j)<YMIN 
                x2_new(i,j)=YMAX-abs(x2_new(i,j)-YMIN);
                4;
              
            end
        end
    end
    x1=x1_new;
    x2=x2_new;
    V1=V1_new;
    V2=V2_new;
    u1=Uo*cos(x1./a).*cos(x2./a); % fluid velocity in x direction
    u2=Uo*sin(x1./a).*sin(x2./a); % fluid velocity in y direction
  
    
    f(t)=figure('visible','off');
      
    for i=1:N
        for j=1:N
        hold on
        scatter(x1(i,j),x2(i,j),5,'k','o'); %Topology
        axis([XMIN-0.5*XMIN XMAX+0.5*XMAX YMIN-0.5*YMIN YMAX+0.5*YMAX]);
        end
    end
    %saveas(f(t),'newout','fig');
    frame(t)=getframe(gcf);
    
end
v=VideoWriter('plot1.avi');
open(v);
for i=1:t
     writeVideo(v,frame(i));
end
close(v);

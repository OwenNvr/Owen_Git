clear all;
clc;

X=[1.5; 2]; %Position at the center of the square

%Guess points
l1=2; %Length link 1
l2=2; %Length link 2
% th1=50*pi/180; %First angle
% th2=50*pi/180; %Second angle
r=0.05; %Radius

x0=[l1;l2;r];

%lower and upper bounds 

lb=[0.1;0.1;0.01];
ub=[10;10;0.5];

%Young module

E=210*(10^9);

%Volumic mass
p=7800;

%Tops of the square
T1=[(X(1)-0.25);(X(2)-0.25)]; %Bottom left
T2=[X(1)+0.25;X(2)-0.25]; %Bottom right
T3=[X(1)+0.25;X(2)+0.25]; %Top right
T4=[X(1)-0.25;X(2)+0.25]; %Top left

%options and fmincon

options = optimset('Display','iter','Tolx',1e-10,'Tolfun',1e-10,'MaxIter',200000,'MaxFunEvals',800000); %Options for the iterations of the gradient descent
[f,val]= fmincon(@(x0)objective5v2(x0,p),x0,[],[],[],[],lb,ub,@(x0)const5v2(x0,E,X),options);

f;
val;

%%%%%%%%%%%For the plot%%%%%%%%%%

th2=(pi-acos(((f(1)^2)+(f(2)^2)-sqrt((X(1)^2)+(X(2)^2))^2)/(2*f(1)*f(2)))); 
th1=atan(X(2)/X(1))-atan(f(2)*sin(th2)/(f(1)+f(2)*cos(th2)));

x1=0+f(1)*cos(th1);
y1=0+f(1)*sin(th1);

x2=x1+f(2)*cos(th1+th2);
y2=y1+f(2)*sin(th1+th2);

plot([0, x1],[0, y1],'r','LineWidth',2);
hold on;
grid on;
plot([x1, x2],[y1, y2],'r','LineWidth',2);

plot([T1(1) T2(1)],[T1(2) T2(2)], 'b','LineWidth',2);
plot([T2(1) T3(1)],[T2(2) T3(2)], 'b','LineWidth',2);
plot([T3(1) T4(1)],[T3(2) T4(2)], 'b','LineWidth',2);
plot([T4(1) T1(1)],[T4(2) T1(2)], 'b','LineWidth',2);

xlim([0,3])
ylim([0,3])
xlabel('x (m)')
ylabel('y (m)')
title('Optimized dimensions of the 2R robot')
txt1 = ['Length l1 : ' num2str(f(1)) ' m'];
txt2 = ['Length l2 : ' num2str(f(2)) ' m'];
txt3 = ['Radius r : ' num2str(f(3)) ' m'];
txt4 = ['Mass m : ' num2str(val) ' kg'];
text(0.2,2,txt1)
text(0.2,1.8,txt2)
text(0.2,1.6,txt3)
text(0.2,1.4,txt4)
hold off;
grid off;


%%%%%Verification of constraints%%%%%%%%%%

 %Jacobian matrix, transpose, and inverse
J=[-f(1)*sin(th1)-f(2)*sin(th1+th2), -f(2)*sin(th1+th2); f(1)*cos(th1)+f(2)*cos(th1+th2), f(2)*cos(th1+th2)];
T_J=transpose(J);
Inv_J=inv(J);

%Product of the Jacobian matrixes
Prod=T_J*J;

%Finding eigen values of J
M=eig(Prod);

%Condition number of the jacobian matrix
K=sqrt(max(M)/min(M));
K


%Calculation of the stiffness matrix
I1=pi*((2*f(3))^4)/64;
I2=pi*((2*f(3))^4)/64;
k11=3*E*I1/(f(1)^3);
k22=3*E*I2/(f(2)^3);

K_th=[k11, 0; 0, k22];

%Calculation of cartesian Stiffness Matrix
K_x=inv(T_J)*K_th*Inv_J;

%Inverse of cartesian matrix multiplied by the force vector
G= inv(K_x)*[100;100];
G






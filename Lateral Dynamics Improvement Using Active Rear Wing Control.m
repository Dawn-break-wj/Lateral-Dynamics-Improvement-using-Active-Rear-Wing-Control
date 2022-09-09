%% 求解jacobian矩阵
%************************************************************************
% clc;
% clear all;
% syms Cf mu m l2 l1 rho Vx Aw Cr Cl Cd a b Iz
% syms Vy omega y phi
% syms delta theta_1 theta_2
% 
% Vy_dot=((Cf+mu*m*(l2/(l1+l2))-mu*((Cl*rho*Vx^2*Aw*theta_1)/2+(Cl*rho*Vx^2*Aw*theta_2)/2)*((b-l2)/(l1+l2)))*(delta-(Vy+l1*omega)/Vx))/m+...
%     ((Cr+mu*m*(l1/(l1+l2))+mu*((Cl*rho*Vx^2*Aw*theta_1)/2+(Cl*rho*Vx^2*Aw*theta_2)/2)*((b+l1)/(l1+l2)))*((-Vy+l2*omega)/Vx))/m-Vx*omega;
% omega_dot=l1*((Cf+mu*m*(l2/(l1+l2))-mu*((Cl*rho*Vx^2*Aw*theta_1)/2+(Cl*rho*Vx^2*Aw*theta_2)/2)*((b-l2)/(l1+l2)))*(delta-(Vy+l1*omega)/Vx))/Iz-...
%     l2*((Cr+mu*m*(l1/(l1+l2))-mu*((Cl*rho*Vx^2*Aw*theta_1)/2+(Cl*rho*Vx^2*Aw*theta_2)/2)*((b+l1)/(l1+l2)))*((-Vy+l2*omega)/Vx))/Iz+...
%     a/Iz*((Cd*rho*Vx^2*Aw*theta_1)/2-(Cd*rho*Vx^2*Aw*theta_2)/2);
% y_dot=Vy;
% phi_dot=omega;
% 
% f=[Vy_dot omega_dot y_dot phi_dot];
% kesi=[Vy omega y phi];
% v=[delta theta_1 theta_2];
% 
% R=jacobian(f,kesi)
% R2=jacobian(f,v) 
% %平衡点处
% Vy=0;omega=0;delta=0;theta_1=0;theta_2=0;
% 
% R1S=subs(R)
% R2S=subs(R2)

%% 对矩阵进行定义和计算
Cf=40000;l1=1.5;l2=1.7;a=0.8;b=1.7;m=1700;mu=0.7;Cr=40000;Iz=2000;Aw=0.24;Cd=1;rho=1.225;
%变量
Vx=55.56;
A=[- (Cf + (l2*m*mu)/(l1 + l2))/(Vx*m) - (Cr + (l1*m*mu)/(l1 + l2))/(Vx*m),    (l2*(Cr + (l1*m*mu)/(l1 + l2)))/(Vx*m) - (l1*(Cf + (l2*m*mu)/(l1 + l2)))/(Vx*m) - Vx, 0, 0;
   (l2*(Cr + (l1*m*mu)/(l1 + l2)))/(Iz*Vx) - (l1*(Cf + (l2*m*mu)/(l1 + l2)))/(Iz*Vx), - (l1^2*(Cf + (l2*m*mu)/(l1 + l2)))/(Iz*Vx) - (l2^2*(Cr + (l1*m*mu)/(l1 + l2)))/(Iz*Vx), 0, 0;
                                                                                   1,                                                                                       0, 0, 0;
                                                                                   0,                                                                                       1, 0, 0];
B=[       (Cf + (l2*m*mu)/(l1 + l2))/m,                         0,                          0;
 (l1*(Cf + (l2*m*mu)/(l1 + l2)))/Iz, (Aw*Cd*Vx^2*a*rho)/(2*Iz), -(Aw*Cd*Vx^2*a*rho)/(2*Iz);
  0, 0, 0;
  0, 0, 0];
Q=[0.01 0 0 0;0 0.01 0 0;0 0 0.01 0;0 0 0 0.01];R=[800 0 0 ;0 1 0;0 0 1];
k=lqr(A,B,Q,R)
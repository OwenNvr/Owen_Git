function [g,geq] = const5v2(x0,E,X)
    %Tops of the square
    T1=[(X(1)-0.25);(X(2)-0.25)]; %Bottom left
    T2=[X(1)+0.25;X(2)-0.25]; %Bottom right
    T3=[X(1)+0.25;X(2)+0.25]; %Top right
    T4=[X(1)-0.25;X(2)+0.25]; %Top left

    %FIRST TOP : BOTTOM LEFT

    th2=(pi-acos(((x0(1)^2)+(x0(2)^2)-sqrt((T1(1)^2)+(T1(2)^2))^2)/(2*x0(1)*x0(2)))); 
    th1=atan(T1(2)/T1(1))-atan(x0(2)*sin(th2)/(x0(1)+x0(2)*cos(th2)));

    %Jacobian matrix, transpose, and inverse
    J=[-x0(1)*sin(th1)-x0(2)*sin(th1+th2), -x0(2)*sin(th1+th2); x0(1)*cos(th1)+x0(2)*cos(th1+th2), x0(2)*cos(th1+th2)];
    T_J=transpose(J);
    Inv_J=inv(J);

    %Product of the Jacobian matrixes
    Prod=T_J*J;

    %Finding eigen values of J
    M=eig(Prod);

    %Condition number of the jacobian matrix
    K=sqrt(max(M)/min(M));

    %CONSTRAINT 1
    g(1)=0.2-K;

    %Calculation of the stiffness matrix
    I1=pi*((2*x0(3))^4)/64;
    I2=pi*((2*x0(3))^4)/64;
    k11=3*E*I1/(x0(1)^3);
    k22=3*E*I2/(x0(2)^3);
    
    K_th=[k11, 0; 0, k22];

    %Calculation of cartesian Stiffness Matrix
    K_x=inv(T_J)*K_th*Inv_J;
    
    %Inverse of cartesian matrix multiplied by the force vector
    G= inv(K_x)*[100;100];
    

    %ONSTRAINT 2
    g(2)=G(1)-0.1*10^(-3);

    %CONSTRAINT 3
    g(3)=G(2)-0.1*10^(-3);

    %SECOND TOP : BOTTOM RIGHT

    th2=(pi-acos(((x0(1)^2)+(x0(2)^2)-sqrt((T2(1)^2)+(T2(2)^2))^2)/(2*x0(1)*x0(2))));
    th1=atan(T2(2)/T2(1))-atan(x0(2)*sin(th2)/(x0(1)+x0(2)*cos(th2)));

    %Jacobian matrix, transpose, and inverse
    J=[-x0(1)*sin(th1)-x0(2)*sin(th1+th2), -x0(2)*sin(th1+th2); x0(1)*cos(th1)+x0(2)*cos(th1+th2), x0(2)*cos(th1+th2)];
    T_J=transpose(J);
    Inv_J=inv(J);

    %Product of the Jacobian matrixes
    Prod=T_J*J;

    %Finding eigen values of J
    M=eig(Prod);

    %Condition number of the jacobian matrix
    K=sqrt(max(M)/min(M));

    %CONSTRAINT 4
    g(4)=0.2-K;

    %Calculation of the stiffness matrix
    I1=pi*((2*x0(3))^4)/64;
    I2=pi*((2*x0(3))^4)/64;
    k11=3*E*I1/(x0(1)^3);
    k22=3*E*I2/(x0(2)^3);
    
    K_th=[k11, 0; 0, k22];

    %Calculation of cartesian Stiffness Matrix
    K_x=inv(T_J)*K_th*Inv_J;
    
    %Inverse of cartesian matrix multiplied by the force vector
    G= inv(K_x)*[100;100];
    

    %ONSTRAINT 5
    g(5)=G(1)-0.1*10^(-3);

    %CONSTRAINT 6
    g(6)=G(2)-0.1*10^(-3);

    %THIRD TOP : TOP RIGHT

    th2=(pi-acos(((x0(1)^2)+(x0(2)^2)-sqrt((T3(1)^2)+(T3(2)^2))^2)/(2*x0(1)*x0(2))));
    th1=atan(T3(2)/T3(1))-atan(x0(2)*sin(th2)/(x0(1)+x0(2)*cos(th2)));

    %Jacobian matrix, transpose, and inverse
    J=[-x0(1)*sin(th1)-x0(2)*sin(th1+th2), -x0(2)*sin(th1+th2); x0(1)*cos(th1)+x0(2)*cos(th1+th2), x0(2)*cos(th1+th2)];
    T_J=transpose(J);
    Inv_J=inv(J);

    %Product of the Jacobian matrixes
    Prod=T_J*J;

    %Finding eigen values of J
    M=eig(Prod);

    %Condition number of the jacobian matrix
    K=sqrt(max(M)/min(M));

    %CONSTRAINT 7
    g(7)=0.2-K;

    %Calculation of the stiffness matrix
    I1=pi*((2*x0(3))^4)/64;
    I2=pi*((2*x0(3))^4)/64;
    k11=3*E*I1/(x0(1)^3);
    k22=3*E*I2/(x0(2)^3);
    
    K_th=[k11, 0; 0, k22];

    %Calculation of cartesian Stiffness Matrix
    K_x=inv(T_J)*K_th*Inv_J;
    
    %Inverse of cartesian matrix multiplied by the force vector
    G= inv(K_x)*[100;100];
    

    %ONSTRAINT 8
    g(8)=G(1)-0.1*10^(-3);

    %CONSTRAINT 9
    g(9)=G(2)-0.1*10^(-3);

    %FOURTH TOP : TOP LEFT

    th2=(pi-acos(((x0(1)^2)+(x0(2)^2)-sqrt((T4(1)^2)+(T4(2)^2))^2)/(2*x0(1)*x0(2))));
    th1=atan(T4(2)/T4(1))-atan(x0(2)*sin(th2)/(x0(1)+x0(2)*cos(th2)));

    %Jacobian matrix, transpose, and inverse
    J=[-x0(1)*sin(th1)-x0(2)*sin(th1+th2), -x0(2)*sin(th1+th2); x0(1)*cos(th1)+x0(2)*cos(th1+th2), x0(2)*cos(th1+th2)];
    T_J=transpose(J);
    Inv_J=inv(J);

    %Product of the Jacobian matrixes
    Prod=T_J*J;

    %Finding eigen values of J
    M=eig(Prod);

    %Condition number of the jacobian matrix
    K=sqrt(max(M)/min(M));

    %CONSTRAINT 10
    g(10)=0.2-K;

    %Calculation of the stiffness matrix
    I1=pi*((2*x0(3))^4)/64;
    I2=pi*((2*x0(3))^4)/64;
    k11=3*E*I1/(x0(1)^3);
    k22=3*E*I2/(x0(2)^3);
    
    K_th=[k11, 0; 0, k22];

    %Calculation of cartesian Stiffness Matrix
    K_x=inv(T_J)*K_th*Inv_J;
    
    %Inverse of cartesian matrix multiplied by the force vector
    G= inv(K_x)*[100;100];
    

    %ONSTRAINT 11
    g(11)=G(1)-0.1*10^(-3);

    %CONSTRAINT 12
    g(12)=G(2)-0.1*10^(-3);

    %Equality constraints
    geq=[];

end





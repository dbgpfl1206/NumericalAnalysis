%memory initialization
clc;

%declaration of variables
m=101; v=1.0;
dx=1.0; dt=0.02; nmax=2500;

%selection of scheme
alpha=input('Input alpha(0, 1.0) ');
theta=input('Input theta(0, 0.5, 1.0) ');
pe=input('Input pe(0.1, 1, 10) ');
scheme=input('Input scheme(1, 2) ');

%grid generation
%when i=1, its distance is at (1*dx)m from origin.
for i=1:m
    x(i)=(i-1)*dx;
end

%initialization
for i=1:m
    rhs(i)=0;
    a1(i)=0;
    b1(i)=0;
    c1(i)=0;
end

%save the initial condition using open file("initial.dat")
fileID=fopen('initial.dat','w');
for i=1:m
    if i==1
        c(i)=1.0;
        fprintf(fileID,'%f, %f',x(i),c(i));
        fprintf(fileID,'\n');
    else
        c(i)=0.0;
        fprintf(fileID,'%f, %f',x(i),c(i));
        fprintf(fileID,'\n');
    end
    c0(i)=c(i); %to plot of initial concentration
end
fclose(fileID);

for iter=1:nmax
    if iter==round(iter/100)*100
        fprintf('%d \n',iter);
    end
    for i=2:m-1 %for2-start statement
        %implicit method
        if scheme==1 %if-start statement
            a1(i)=(-v*alpha*theta/(2*dx)-v*(1-alpha)*theta/dx...
                -v*theta/(pe*dx));
            b1(i)=(1/dt+v*(1-alpha)*theta/dx+2*v*theta/(pe*dx));
            c1(i)=(v*alpha*theta/(2*dx)-v*theta/(pe*dx));
            d1=(v*alpha*(1-theta)/(2*dx)+v*(1-alpha)*(1-theta)/dx...
                +v*(1-theta)/(pe*dx));
            e1=(1/dt-v*(1-alpha)*(1-theta)/dx-2*v*(1-theta)/(pe*dx));
            f1=(-v*alpha*(1-theta)/(2*dx)+v*(1-theta)/(pe*dx));
            rhs(i)=d1*c(i-1)+e1*c(i)+f1*c(i+1);
        %explicit method
        elseif scheme==2
            d1=(v*alpha*(1-theta)/(2*dx)+v*(1-alpha)*(1-theta)/dx...
                +v*(1-theta)/(pe*dx));
            e1=(1/dt-v*(1-alpha)*(1-theta)/dx-2*v*(1-theta)/(pe*dx));
            f1=(-v*alpha*(1-theta)/(2*dx)+v*(1-theta)/(pe*dx));
            c(i)=d1*c(i-1)+e1*c(i)+f1*c(i+1);
            c(i)=c(i)*dt;
        end %if-end statement
    end %for2-end statement
    
    %implicit method
    if scheme==1 %if-implicit start statement
        rhs(2)=rhs(2)-a1(2)*c(1);
        rhs(m-1)=rhs(m-1)-c1(m-1)*c(m);
        %solve the Thomas algorithm
        for i=3:m-1
            r=a1(i)/b1(i-1);
            b1(i)=b1(i)-r*c1(i-1);
            rhs(i)=rhs(i)-r*rhs(i-1);
        end
        rhs(m-1)=rhs(m-1)/b1(m-1);
        for i=m-2:-1:2
            rhs(i)=(rhs(i)-c1(i)*rhs(i+1))/b1(i);
        end %end of Thomas algorithm
        
        %renewal of concentration at the next time step
        for i=2:m-1
            c(i)=rhs(i);
        end
    end %if-implicit end statement
    
    %renewal of boundary conditions
    c(1)=1;       
    c(m)=c(m-1);
end %for1-start statement

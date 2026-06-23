function dX = int_dyn(X,label,time,sig,t,flag)

if strcmp(flag,'Train')==1
u = sig;
end


if max(abs(sig))>0
if strcmp(flag,'Simulate')
u  = interp1(time,sig,t,'nearest');
end
else
    u = 0;
end


%% 1D Systems 
if  strcmp(label,'Integrator')==1
dX(1,:) = u;
end

%% 3D Chaotic Attractors 
if strcmp(label,'Rossler')==1 
a = 0.2;
b = 0.2;
c = 5.7;
TX = 10;

x = X(1,:);
y = X(2,:); 
z = X(3,:);

dX(1,:)= -y-z;
dX(2,:)= x+a*y;
dX(3,:)= b+z.*(x-c);
dX  = dX/TX;
end

if strcmp(label,'Lorenz')==1
rho = 28;
sigma = 10;
beta = 8/3;
TX = 40;

x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=sigma*(y-x);
dX(2,:)=x.*(rho-z)-y;
dX(3,:)=x.*y-beta*z;
dX = dX/TX;
end

if strcmp(label,'Thomas')
    TX = 1;
b = 0.208186;    
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=sin(y)-b*x;
dX(2,:)=sin(z)-b*y;
dX(3,:)=sin(x)-b*z;
dX = dX/TX;
end


if strcmp(label,'SprottA')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y;
dX(2,:)=-x+y.*z;
dX(3,:)=1-y.^2;
dX = dX/TX;
end


if strcmp(label,'SprottB')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y.*z;
dX(2,:)=x-y;
dX(3,:)=1-x.*y;
dX = dX/TX;
end


if strcmp(label,'SprottC')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y.*z;
dX(2,:)=x-y;
dX(3,:)=1-x.^2;
dX = dX/TX;
end

if strcmp(label,'SprottD')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-y;
dX(2,:)=x+z;
dX(3,:)=x.*z+3*y.^2;
dX = dX/TX;
end


if strcmp(label,'SprottE')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y.*z;
dX(2,:)=x.^2-y;
dX(3,:)=1-4*x;
dX = dX/TX;
end


if strcmp(label,'SprottF')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y+z;
dX(2,:)=-x+0.5*y;
dX(3,:)=x.^2-z;
dX = dX/TX;
end

if strcmp(label,'SprottG')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=0.4*x+z;
dX(2,:)=x.*z-y;
dX(3,:)=-x+y;
dX = dX/TX;
end

if strcmp(label,'SprottH')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-y+z.^2;
dX(2,:)=x+0.5*y;
dX(3,:)=x-z;
dX = dX/TX;
end


if strcmp(label,'SprottI')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-0.2*y;
dX(2,:)=x+z;
dX(3,:)=x+y.^2-y;
dX = dX/TX;
end

if strcmp(label,'SprottI')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-0.2*y;
dX(2,:)=x+z;
dX(3,:)=x+y.^2-z;
dX = dX/TX;
end

if strcmp(label,'SprottJ')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=2*z;
dX(2,:)=-2*y+z;
dX(3,:)=-x+y+y.^2;
dX = dX/TX;
end

if strcmp(label,'SprottK')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=x.*y - z;
dX(2,:)=x-y;
dX(3,:)=x+0.3*z;
dX = dX/TX;
end


if strcmp(label,'SprottL')
TX = 1;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y+3.9*z;
dX(2,:)=0.9*(x.^2)-y;
dX(3,:)=1-x;

dX = dX/TX;
end


if strcmp(label,'SprottM')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-z;
dX(2,:)=-x.^2 - y;
dX(3,:)=1.7+1.7*x  + y;
dX = dX/TX;
end

if strcmp(label,'SprottN')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-2*y;
dX(2,:)=x+z.^2;
dX(3,:)= 1+y-2*z;
dX = dX/TX;
end

if strcmp(label,'SprottO')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=y;
dX(2,:)=x-z;
dX(3,:)= x+x.*z+2.7*y;
dX = dX/TX;
end

if strcmp(label,'SprottP')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=2.7*y + z;
dX(2,:)=-x+y.^2;
dX(3,:)= x+y;
dX = dX/TX;
end

if strcmp(label,'SprottQ')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-z;
dX(2,:)=x-y;
dX(3,:)= 3.1*x+y.^2 + 0.5*z;
dX = dX/TX;
end

if strcmp(label,'SprottR')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=0.9-y;
dX(2,:)=0.4+z;
dX(3,:)= x.*y-z;
dX = dX/TX;
end

if strcmp(label,'SprottS')
TX = 10;   
x=X(1,:);
y=X(2,:);
z=X(3,:);

dX(1,:)=-x-4*y;
dX(2,:)=x+z.^2;
dX(3,:)=1+x;
dX = dX/TX;
end



if strcmp(label,'Vanderpol')
TX = 10;  
mu = 5;
x=X(1,:);
y=X(2,:);

dX(1,:)=mu*(x-(x.^3)/3 - y);
dX(2,:)=x/mu;

dX = dX/TX;
end



if strcmp(label,'Pitchfork')
TX = 1;  
x=X(1,:);


dX(1,:)=0.5*x-x.^3;
dX = dX/TX;
end

if strcmp(label,'Hopf Normal Form')
beta = 0.5;
sigma =  -1;
x=X(1,:);
y=X(2,:);
dX(1,:) = beta*x  - y + sigma*x.*(x.^2 + y.^2);
dX(2,:) = x + beta*y + sigma*y.*(x.^2 + y.^2); 
end


if strcmp(label,'Chua1')
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = 0.3*y+x-x.^3;
dX(2,:) = x+z;
dX(3,:) = -y;
dX = dX/TX;
end


if strcmp(label,'Chua2')
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = 0.2*y + -x +  2*tanh(x);
dX(2,:) = x+z;
dX(3,:) = -y;
dX = dX/TX;
end

if strcmp(label,'Chua3')
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = 0.2*y +x -x.*abs(x);
dX(2,:) = x+z;
dX(3,:) = -y;
dX = dX/TX;
end

if strcmp(label,'Chua4')
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = 0.2*y -x -2*sin(x);
dX(2,:) = x+z;
dX(3,:) = -y;
dX = dX/TX;
end

if strcmp(label,'Chua5')
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = 0.2*y -0.3*x + sign(x);
dX(2,:) = x+z;
dX(3,:) = -y;
dX = dX/TX;
end

if strcmp(label,'Chua6')
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = 0.2*y -  x  + 2*atan(x);
dX(2,:) = x+z;
dX(3,:) = -y;
dX = dX/TX;
end




if strcmp(label,'Rikitake')
mu  = 1;
alpha = 1;
TX  = 10;
x =  X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = -mu*x + y.*z;
dX(2,:) = -mu*y + x.*(z-alpha);
dX(3,:) = 1-x.*y;
dX = dX/TX;
end


if strcmp(label,'Nose Hoover')
TX  = 10;
x = X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = y;
dX(2,:) = y.*z - x;
dX(3,:) = 1-y.^2;
dX = dX/TX;
end



if strcmp(label,'Halvorsen')
TX  = 10;
x = X(1,:);
y = X(2,:);
z = X(3,:); 

dX(1,:) = -a*x - 4*y - 4*z - y.^2;
dX(2,:) = -a*y - 4*z - 4*x - z.^2;
dX(3,:) =  a*z - 4*x - 4*y - x.^2;
dX = dX/TX;
end


if strcmp(label,'MO0')
TX  = 10;
a = 0.6;
b= 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = abs(x)-1;
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO1')
TX  = 10;
a = 0.6;
b= 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = 1 - 6*max(x,0);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO2')
TX  = 10;
a = 0.6;
b= 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = sign(x)-x;
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO3')
TX  = 10;
a = 1;
b= 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = 1.1*(x.^2-1);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO4')
TX  = 10;
a = 0.5;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = x.*(x-1);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO5')
TX  = 10;
a = 0.7;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = x.*(1-x.^2);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO6')
TX  = 10;
a = 0.4;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = (x.^2).*(1-x);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO7')
TX  = 10;
a = 0.6;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = (x.^2).*(1-x.^2);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO8')
TX  = 10;
a = 0.5;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = (x).*(x.^4-1);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO9')
TX  = 10;
a = 0.4;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = (x.^3).*(1-x);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO10')
TX  = 10;
a = 0.6;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = (x.^2).*(1-x.^3);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO11')
TX  = 10;
a = 1;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = 5 - exp(x);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO12')
TX  = 10;
a = 1;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g =  7 - 8*tanh(x);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO13')
TX  = 10;
a = 1;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = 6*tanh(x) - 3*x;
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end

if strcmp(label,'MO14')
TX  = 10;
a = 0.6;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = 6*atan(x) - x;
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end


if strcmp(label,'MO15')
TX  = 10;
a = 0.6;
b = 1;
x = X(1,:);
y = X(2,:);
z = X(3,:); 


g = x - 0.5*sinh(x);
dX(1,:) = y; 
dX(2,:) = z;
dX(3,:) = g - a*z - b*y;
dX = dX/TX;
end




















end
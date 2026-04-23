clear all
close all
clc
tic
%%
net_config %confiure initial network. 
epsilon = 0.2e-1; %gradient descent rate. 
no = 3*10^5; %number of iterates in gradient descent. 
%% Run a simulation of the unscaled dynamical system to determine rescaling onto the unit cube. 
load network_configuration.mat
k0 = size(phi,2);
for learn_sys =1:40

%% generate training points. 
x = 2*rand(k0,nx)-1;   %random sample points across space for training. 
xv =  2*rand(k0,100)-1; %random sample points across space for validation.
c = 2*rand(1,nx)-1;
cv = 2*rand(1,100)-1;

dt = 1e-2;    
k = dim(learn_sys);
bias0 = randn(N,1);% reinitialize the bias current.     
y0 =  ic(learn_sys,1:k)';
Ts = 1000; %sample time to determine how to scale dynamics onto unit cube.
ns = round(Ts/dt);
times = (0:ns)*dt;
sig = cos(2*pi*times'*40*rand(1,20)/1000)*rand(20,1);
sig(abs(sig)>1)=0;

flag = 'Simulate';
tic
label = dynamics{learn_sys};
[t,y] = ode45(@(t,y) int_dyn(y,label,times,input(learn_sys,1)*sig,input(learn_sys,1)*t,flag),[0,Ts],y0);  %integrate. 
disp(sprintf('Running Sample of %s',label))
toc
xk = x(1:k,:);
xvk = xv(1:k,:);

%%
flag = 'Train';
if rescale(learn_sys)==1
scal_mat = diag(max(abs(y))); %estimate rescaling factor.
mean_mat = mean(y)'; %estimate means. 
xscaled = scal_mat*xk + mean_mat; %change to mean  0 and rescale to max 1 
xscaledv = scal_mat*xvk + mean_mat; %change to mean  0 and rescale to max 1
gx = pinv(scal_mat)*int_dyn(xscaled,label,times,input(learn_sys,1)*c,t,flag)+ xk;  %determine the scaled derivative. 
gxv = pinv(scal_mat)*int_dyn(xscaledv,label,times,input(learn_sys,1)*cv,t,flag) + xvk;  %determine the scaled derivative. 
else
scal_mat = eye(k); 
mean_mat =  0;
xscaled = scal_mat*xk + mean_mat; %change to mean  0 and rescale to max 1 
xscaledv = scal_mat*xvk + mean_mat; %change to mean  0 and rescale to max 1
gx = pinv(scal_mat)*int_dyn(xscaled,label,times,input(learn_sys,1)*c,t,flag)+ xk;  %determine the scaled derivative. 
gxv =pinv(scal_mat)*int_dyn(xscaledv,label,times,input(learn_sys,1)*cv,t,flag)+ xvk;  %determine the scaled derivative. 
end



%% Run a simulation of the scaled system for comparison at the end.
flag = 'Simulate';
T = 1000;
nt = round(T/dt);
timet = (0:nt)*dt;
sigt = cos(2*pi*timet'*40*rand(1,20)/1000)*rand(20,1)/500;
tic
[t,y] = ode45(@(t,y) pinv(scal_mat)*int_dyn(scal_mat*y+mean_mat,label,timet,input(learn_sys,1)*sigt,t,flag),0:dt:T,pinv(scal_mat)*(y0 - mean_mat)); % system with the appropriately scaled dynamics. 
disp(sprintf('Running Rescaled/Comparison %s',label))
toc
%% Run gradient descent (or variants).  
 if k<3  
 i1 = find(sum(abs(eta(:,k+1:end)),2)==0); %find the neurons that only encode in the first k dimensions. 
 i2 = find(sum(abs(eta(:,k+1:end)),2)>0); %index of all the other neurons. 
 else 
  i1 = 1:N;
  i2 = [];
 end

%%
tic
[bias,xhat,store,storeb] = bias_gradient_descent(xk,xvk,input(learn_sys,1)*c,input(learn_sys,1)*cv,win,eta,bias0,phi,gx,gxv,no,epsilon,i1,i2);
disp(sprintf('Learning Bias Currents for %s',label))
toc
f1 = figure(1)
set(f1,'position',[0,0,2000,1000])
print(f1,sprintf('gradient_%s.jpg',dynamics{learn_sys}),'-djpeg','-r300')
close(f1);
%% Integrate the network.
y0net = [pinv(scal_mat)*(y0 - mean_mat);0.1*randn(size(phi,2)-k,1)];
[t1,y1] = ode45(@(t,y) net_int(eta,phi,bias,win,input(learn_sys,1)*sigt,timet,t,y),0:dt:T,y0net);
f2 = figure(learn_sys*10);
plot(t,y,'r'),  hold on
plot(t1,y1(:,1:k),'k--')
title(dynamics{learn_sys})
f2 = figure(learn_sys*10);
set(f2,'position',[0,0,2000,1000])
print(f2,sprintf('net_%s.jpg',dynamics{learn_sys}),'-djpeg','-r300')
f3 = figure(learn_sys*100);
subplot(1,2,1)
plot(y1(:,1),y1(:,2))
xlim([-1.1,1.1])
ylim([-1.1,1.1])
title('Network')
if k>=2
subplot(1,2,2)
plot(y(:,1),y(:,2))
xlim([-1.1,1.1])
ylim([-1.1,1.1])
title(dynamics{learn_sys})
set(f3,'position',[0,0,2000,1000])
print(f3,sprintf('phase_%s.jpg',dynamics{learn_sys}),'-djpeg','-r300')
end


save(sprintf('%s_results.mat',dynamics{learn_sys}),'bias','t1','y1','eta','phi','t','y','store','storeb','mean_mat','scal_mat','epsilon','no')
close all
end

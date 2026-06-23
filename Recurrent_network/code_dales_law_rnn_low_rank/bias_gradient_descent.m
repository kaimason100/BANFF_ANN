function [bias,xhat,store,storeb] = bias_gradient_descent(x,xv,c,cv,win,eta,bias,phi,gx,gxv,no,epsilon,i1,i2) 
k = size(gx,1);
phik = size(phi,2);
store = zeros(no,2);
storeb = zeros(no,10);


x(k+1:phik,:) = 0;
xv(k+1:phik,:) = 0;

prod1 = eta*x;
prod2 = eta*xv;
theta =  zeros(length(i1),1);
gamma = 0.9;


prod3 = win*c;
prod4 = win*cv;
bias(i2) = -20;


for j = 1:no 
fx = tanh(prod1+bias+prod3);
fxv = tanh(prod2+bias+prod4);
xhatv = phi(:,1:k)'*(fxv) ;
xhat = phi(:,1:k)'*(fx) ; 
dnu = (1-fx(i1,:).^2);

err =  (xhat - gx)';
errv = (xhatv-gxv)';

grad = sum(phi(i1,1:k).*(dnu*err),2);
theta =  gamma*theta + epsilon*grad;
bias(i1,1) = bias(i1,1)  - theta; 
loss = mean(mean(err.^2));
lossv =  mean(mean(errv.^2));

if loss>10^4
break;
end
store(j,:)=[loss,lossv];
storeb(j,:) = bias(i1(1:10));
figure(1)
if mod(j,10)==1
if size(x,1) >=2 
 subplot(1,3,1)
loglog((1:j),(store(1:j,:)))
title('Loss')
xlabel('Iterate')
subplot(1,3,2)
plot(xhat',gx','.'),  hold on 
plot([-1,1],[-1,1],'b'), hold off
xlabel('\hat{y}')
xlabel('y')
title(sprintf('Step %d',j))
subplot(1,3,3)
plot((1:j),storeb(1:j,:))
title('Biases vs. Iterate')
end


end
drawnow

end
end




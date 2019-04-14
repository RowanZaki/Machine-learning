clc
clear all
close all

[num,text] = xlsread('house_prices_data_training_data.csv');
x = num(:,4:end);
m = 18;

% Principal component analysis %
% Step 1 %
Corr_x = corr(x);
% Step 3 %
x_cov=cov(x) ;
% Step 4 %
[U S V] =  svd(x_cov);
% Step 5 %

k = 0;
alpha = 10^5;
while (alpha>0.001) && k<m 
     k = k+1;
     num_k = max(S(:,1:k));
     denom_m = max(S);
     alpha = 1-(sum(num_k)/sum(denom_m));
 
end
%Step 6
Reduced_Data = U(:,1:k)'*x';
%Step 7
X_approximate = U(:,1:k)*Reduced_Data;


% scaling all features of x %
for i=1:length(x(1,:))
    x(:,i) = (x(:,i)-mean(x(:,i)))/std(x(:,i)); 
end
x_approximate=x;

% scaling all features 0f the approximate Xapprox %
for i=1:length(x(1,:))
    X_approximate(:,i) = (X_approximate(:,i)-mean(X_approximate(:,i)))/std(X_approximate(:,i)); 
end

% Step 8
error = (1/m)*(sum(x - X_approximate'));

% Linear regression %

alpha=0.03;
MSE=[];

% normalizing the price %
y = (num(:,3)-mean(num(:,3)))/std(num(:,3)); 
m = length(y);
% generation of x0 %
Variables = [ones(m,1) x_approximate];
n = length(Variables(1,:));
% generation of thetas %
theta_old= ones(1,n);

% generation of h and cost function %
h=[];
for i=1:m
h(i) = theta_old*Variables(i,:)';
end
J_old = (1/(2*m))*sum((h' -y).^2);
% generation of new thetas %
for i=1:n
theta_new= theta_old - (alpha*(1/m)*sum( (h' - y).*Variables(:,i)));
end

for i=1:m
h(i) = theta_new*Variables(i,:)';
end

J_new = (1/(2*m))*sum((h' -y).^2);
iterations =1;
 mse1 = J_new ;
 hu=[];
while J_old - J_new > 10^-1

theta_old = theta_new;
%generation of h and cost function

for i=1:m
h(i) = theta_old*Variables(i,:)';
end
J_old = (1/(2*m))*sum((h' -y).^2);
for i=1:n
theta_new= theta_old - (alpha*(1/m)*sum( (h' - y).*Variables(:,i)));
end
for i=1:m
hu(i) = theta_new*Variables(i,:)';
end
J_new = (1/(2*m))*sum((hu' -y).^2);
mse(iterations) = J_new;
iterations= iterations+1
end
mse = [mse1 mse];
iterations =length(mse);
figure(1)
plot([1:iterations],mse)
xlabel('No. of iterations')
ylabel('MSE')
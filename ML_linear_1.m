clc
clear all
close all

[num,text] = xlsread('house_prices_data_training_data.csv');
x = num(:,4:end);
% scaling all features %
for i=1:length(x(1,:))
    x(:,i) = (x(:,i)-mean(x(:,i)))/std(x(:,i)); 
end

% dividing the data %
x_train = x(1:10800,4:end);
x_cv = x(10800:14400,4:end);
x_test = x(14400:end,4:end);

m_train = length(num(1:10800,1));
m_cv = length(num(10800:14400,1));
m_test = length(num(14400:end,1));

% TRAINING SET %

% creating x0 %
x0 = ones(length(num(1:10800,1)),1);

% linear hypothesis %
% x1_train = x(1:10800,4:11);

% HYPOTHESIS 2 %
x1_train = [x(1:10800,4:8) x(1:10800,4:8).^2];

% HYPOTHESIS 3 %
% x1_train = [x(1:10800,4:12) x(1:10800,4:12).^2 x(1:10800,4:12).^3] ;

% HYPOTHESIS 4 %
% x1_train = [x(1:10800,[4,5,8,11]) x(1:10800,[4,5,8,11]).^2 x(1:10800,[4,5,8,11]).^3 x(1:10800,[4,5,8,11]).^4];


% To normalize y which is the price %
y_train = (num(1:10800,3)-mean(num(1:10800,3)))/std(num(1:10800,3));

% thetas for hyp 1 %
% Theta_old1_train = rand(1,9);

% thetas for hyp 2 %
Theta_old1_train = rand(1,11);

% thetas for hyp 3%
% Theta_old1_train = rand(1,25);

% thetas for hyp 4 %
% Theta_old1_train = rand(1,17);

alpha = 0.003;
num_itr = 100;

%Hyp 1%
% xnorm_train = [x0 x1_train];
% n_train = length(xnorm_train(1,:));

%Hyp 2%
xnorm_train = [x0 x1_train];
% n_train = length(xnorm_train(1,:));

%Hyp 3%
% xnorm_train = [x0 x1_train];
% n_train = length(xnorm_train(1,:));

%Hyp 4%
% xnorm_train = [x0 x1_train];
% n_train = length(xnorm_train(1,:));


% TRAINING SET %
diff_train = [];
for j=1:length(xnorm_train(:,1))
    h_train = sum(Theta_old1_train.*xnorm_train(j,:));
    diff_train = [diff_train;h_train-y_train(j)];
end

mse_tot=[];
% To update thetas %

% for each hypothesis i loop on the four different values of alpha to see which one is better
% for this hypo %
for k=1:length(alpha)
    % to empty the mse for each alpha %
    mse = [];
for l=1:num_itr
  % 4 new thetas %  
% Theta_new1_train = zeros(1,9);

% Hypo  2%
Theta_new1_train = zeros(1,11);

% Hypo 3 %
% Theta_new1_train = zeros(1,25);

% Hypo 4 %
% Theta_new1_train = zeros(1,17);

for j=1:length(Theta_old1_train)
    Theta_new1_train(1,j) = Theta_old1_train(1,j)- alpha(k)*(1/m_train)*(sum(diff_train.*xnorm_train(:,j)));
end
Theta_old1_train = Theta_new1_train;
diff_train = [];
for j=1:length(xnorm_train(:,1))
    h_train = sum(Theta_old1_train.*xnorm_train(j,:));
    diff_train = [diff_train;h_train-y_train(j)];
end
mse = [mse MSE(diff_train)];
end
mse_tot = [mse_tot;mse];
end

% cross validation %

y_cv = (num(10800:14400,3)-mean(num(10800:14400,3)))/std(num(10800:14400,3));

% creating x0 %
x0_cv = ones(length(num(10800:14400,1)),1);

% 1 %
% x1_cv = x(10800:14400,4:11);
% xnorm_cv = [x0_cv x1_cv];

% 2 %
x1_cv = [x(10800:14400,4:8) x(10800:14400,4:8).^2];
xnorm_cv = [x0_cv x1_cv];

% 3 %
%x1_cv = [x(10800:14400,4:12) x(10800:14400,4:12).^2 x(10800:14400,4:12).^3] ;
% xnorm_cv = [x0_cv x1_cv];

% 4 %
%x1_cv = [x(10800:14400,[4,5,8,11]) x(10800:14400,[4,5,8,11]).^2 x(10800:14400,[4,5,8,11]).^3 x(10800:14400,[4,5,8,11]).^4];
% xnorm_cv = [x0_cv x1_cv];

diff_cv = [];
for j=1:length(x1_cv(:,1))
    h_cv = sum(Theta_old1_train.*xnorm_cv(j,:));
    diff_cv = [diff_cv;h_cv-y_cv(j)];
    end
mse_cv =MSE_CV(diff_cv);

plot([1:num_itr(1)],mse)
hold on 
plot (mse_cv)

% Normal equation %
%Theta = Normal_equation(xnorm,y);



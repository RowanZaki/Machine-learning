clc
clear all
close all

[num,text] = xlsread('heart_DD.csv');

% To normalize x%
x = num(:,1:13);
for i=1:length(x(1,:))
    x(:,i) = (x(:,i)-mean(x(:,i)))/std(x(:,i)); 
end

y_train = num(1:150,14);
y_cv = num(150:200,14);
y_test = num(200:end,14);


m_train = length(num(1:150,1));
m_cv = length(num(150:200,1));
m_test = length(num(200:end,1));

% HYPO 1 %
x1 = x(1:150,1:8);
n = length(x1(1,:));

% HYPO 2 %
%x1 = [x(1:150,1:4) x(1:150,1:4).^2];
% n = length(x1(1,:));

% HYPO 3 %
%x1 = [x(1:150,6:12) x(1:150,6:12).^2 x(1:150,6:12).^3] ;
% n = length(x1(1,:));

% HYPO 4 %
%x1 = [x(1:150,[4,5,8,11]) x(1:150,[4,5,8,11]).^2 x(1:150,[4,5,8,11]).^3 x(1:150,[4,5,8,11]).^4];
% n = length(x1(1,:));

% creating x0 %
x0 = ones(length(num(1:150,1)),1);


xnorm_train = [x0 x1];

% xnorm_train = [x0 x1];
% xnorm_train = [x0 x1];
% xnorm_train = [x0 x1];


% thetas for hyp 1 %
Theta1_train = zeros(n+1,1);

% thetas for hyp 2 %
% Theta1_train = zeros(n+1,1);

% thetas for hyp 3%
%Theta1_train = zeros(n+1,1);

% thetas for hyp 4 %
%Theta1_train = zeros(n+1,1);


alpha = 0.01;
lamda = 10;
E_train = [];
u = 1;
for i=1:length(Theta1_train)
e = exp(-xnorm_train*Theta1_train);
h = 1./(1+(e));
g(i) = (1/m_train)*sum((h-y_train)'*xnorm_train(:,i));
u = u+1;
Theta1_train = Theta1_train-(alpha/m_train)*xnorm_train'*(h-y_train);
E_train(u) = -(1/m_train)*sum((y_train.*(log(h)))+(1-y_train).*(log(1-h)))+(lamda/(2*m_train))*sum(Theta1_train.^2);

end
u = 1;


% CROSS VALIDATION %

% HYPO 1 %
x1_cv = x(150:200,1:8);
n_cv = length(x1_cv(1,:));

% HYPO 2 %
%x1_cv = [x(150:200,1:4) x(150:200,1:4).^2];
% n_cv = length(x1_cv(1,:));

% HYPO 3 %
%x1_cv = [x(150:200,6:12) x(150:200,6:12).^2 x(150:200,6:12).^3] ;
% n_cv = length(x1_cv(1,:));

% HYPO 4 %
%x1_cv = [x(150:200,[4,5,8,11]) x(150:200,[4,5,8,11]).^2 x(150:200,[4,5,8,11]).^3 x(150:200,[4,5,8,11]).^4];
% n_cv = length(x1_cv(1,:));

% creating x0 %
x0_cv = ones(length(num(150:200,1)),1);


xnorm_cv = [x0_cv x1_cv];

% xnorm_cv = [x0_cv x1_cv];
% xnorm_cv = [x0_cv x1_cv];
% xnorm_cv = [x0_cv x1_cv];


% thetas for hyp 1 %
Theta1_cv = zeros(n_cv+1,1);

% thetas for hyp 2 %
% Theta1_cv = zeros(n_cv+1,1);

% thetas for hyp 3%
%Theta1_cv = zeros(n_cv+1,1);

% thetas for hyp 4 %
%Theta1_cv = zeros(n_cv+1,1);

E_cv = [];
for i=1:length(Theta1_cv)
e_cv = exp(-xnorm_cv*Theta1_cv);
h_cv = 1./(1+(e_cv));
g_cv(i) = (1/m_cv)*sum((h_cv-y_cv)'*xnorm_cv(:,i));
u = u+1;
Theta1_cv = Theta1_cv-(alpha/m_cv)*xnorm_cv'*(h_cv-y_cv);
E_cv(u) = -(1/m_cv)*sum((y_cv.*(log(h_cv)))+(1-y_cv).*(log(1-h_cv)))+(lamda/(2*m_cv))*sum(Theta1_cv.^2);

end
u = 1;


% TESTING %

% HYPO 1 %
x1_test = x(200:end,1:8);
n_test = length(x1_test(1,:));

% HYPO 2 %
%x1_test = [x(200:end,1:4) x(200:end,1:4).^2];
% n_test = length(x1_test(1,:));

% HYPO 3 %
%x1_cv = [x(150:200,6:12) x(150:200,6:12).^2 x(150:200,6:12).^3] ;
% n_cv = length(x1(1,:));

% HYPO 4 %
%x1_test = [x(200:end,[4,5,8,11]) x(200:end,[4,5,8,11]).^2 x(200:end,[4,5,8,11]).^3 x(200:end,[4,5,8,11]).^4];
% n_test = length(x1_test(1,:));

% creating x0 %
x0_test = ones(length(num(200:end,1)),1);


xnorm_test = [x0_test x1_test];

% xnorm_test = [x0_test x1_test];
% xnorm_test = [x0_test x1_test];
% xnorm_test = [x0_test x1_test];


% thetas for hyp 1 %
Theta1_test = zeros(n_test+1,1);

% thetas for hyp 2 %
% Theta1_test = zeros(n_test+1,1);

% thetas for hyp 3%
%Theta1_test = zeros(n_test+1,1);

% thetas for hyp 4 %
%Theta1_test = zeros(n_test+1,1);

E_test = [];
for i=1:length(Theta1_test)
e_test = exp(-xnorm_test*Theta1_test);
h_test = 1./(1+(e_test));
g_test(i) = (1/m_test)*sum((h_test-y_test)'*xnorm_test(:,i));
u = u+1;
Theta1_test = Theta1_test-(alpha/m_test)*xnorm_test'*(h_test-y_test);
E_test(u) = -(1/m_test)*sum((y_test.*(log(h_test)))+(1-y_test).*(log(1-h_test)))+(lamda/(2*m_test))*sum(Theta1_test.^2);

end
u = 1;


plot(E_train)









% for l=1:num_itr
%   % 4 new thetas %  
% Theta_new1_train = zeros(n+1,1);
% %Theta_new1_train = zeros(n+1,1);
% %Theta_new1_train = zeros(n+1,1);
% %Theta_new1_train = zeros(n+1,1);
% 
% for j=1:length(Theta_old1_train)
%     Theta_new1_train(1,j) = Theta_old1_train(1,j)- alpha*(1/m_train)*(sum(diff_train.*xnorm_train(:,j)));
% end
% Theta_old1_train = Theta_new1_train;
% diff_train = [];
% for j=1:length(xnorm_train(:,1))
%     h = 1./(1+exp(-x1*Theta_old1_train));
%     diff_train = [diff_train;h_train-y_train(j)];
% end
% mse = [mse MSE(diff_train)];
% end
% mse_tot = [mse_tot;mse];
% 
% 
% % cross validation %
% 
% % m_cv = length(num(10800:14400,1));
% % y_cv = (num(10800:14400,3)-mean(num(10800:14400,3)))/std(num(10800:14400,3));
% % 
% % % creating x0 %
% % x0_cv = ones(length(num(10800:14400,1)),1);
% % % linear hypothesis %
% % x1_cv = x(10800:14400,4:11);
% % xnorm_cv = [x0_cv x1_cv];
% % 
% % diff_cv = [];
% % for j=1:length(x1_cv(:,1))
% %     h_cv = sum(Theta_old1_train.*xnorm_cv(j,:));
% %     diff_cv = [diff_cv;h_cv-y_cv(j)];
% %     end
% % mse_cv =MSE_CV(diff_cv);






clear all
close all
clc
rng('shuffle'); % Shuffle the random number generator seed

dim = 41;
dim_T=40;

n = 2000;   %NMF
eig_rho =0.682;
W_in_a =0.0235;
alpha = 0.572;
beta = 10^(-3.9404);
k = 0.2535;
bias=1.4794;





%% data simulation

  load NMF_1.mat
  data=data(:,1:1100);
  %data=normalize(data);

% set the washup, training, and testing length
transient_length = 0e3;
washup_length = 0.1e3;
train_length = 0.55e3 + washup_length;
valtesting_length = 50;
testing_length = 300;
short_prediction_length = 50;


tmax = (transient_length + washup_length + train_length + valtesting_length+testing_length + 1);
data_length_all = washup_length + train_length +valtesting_length+ testing_length + 1;


%% reservoir computer configuration

W_in = W_in_a*(2*rand(n,dim)-1);
% hidden network: symmeric, normally distributed, with mean 0 and variance 1.
A = sprandsym(n, k); 
%rescale eig
eig_D=eigs(A,1); %only use the biggest one. Warning about the others is harmless
A=(eig_rho/(abs(eig_D))).*A;
A=full(A);

%% train
disp('training...')
tic;

udata = data;

r_train = zeros(n,train_length - washup_length);
y_train = zeros(dim_T,train_length - washup_length);
r_end = zeros(n,1);

udata=udata';
% train_x: input; train_y: target
train_x = zeros(train_length,dim);
train_y = zeros(train_length,dim_T);
train_x(:,:) = udata(1:train_length,:);
train_y(:,:) = udata(2:train_length+1,1:dim_T);
train_x = train_x';
train_y = train_y';
% initialize r state with all zeros, also can be chosen as random numebrs
r_all = zeros(n,train_length+1);
% update the r state step by step
for ti=1:train_length
    r_all(:,ti+1) = (1-alpha)*r_all(:,ti) + alpha*tanh(A * r_all(:,ti)+W_in*train_x(:,ti)+bias);
end

r_out=r_all(:,washup_length+2:end); 
r_out(2:2:end,:)=r_out(2:2:end,:).^2;
r_end(:)=r_all(:,end); 

r_train(:,:) = r_out;
y_train = train_y(:,washup_length+1:end);
% calculate the weights of Wout
Wout = y_train *r_train'*(r_train*r_train'+beta*eye(n))^(-1);
pred_train=Wout*r_train;

t_train = toc;
disp(['training time: ', num2str(t_train)]);

%% test
disp('testing...')

valtesting_start = train_length+2;

test_pred_y = zeros(valtesting_length,dim_T);
test_real_y = zeros(valtesting_length,1);
test_real_y = udata(valtesting_start:(valtesting_start+valtesting_length-1),1:dim_T);

r=r_end;
% u = zeros(1,1);
%u = udata(train_length+1,1);
u = zeros(1,dim);
u = udata(train_length-1,1:dim);
u=u';

for t_i = 1:valtesting_length

    r = (1-alpha) * r + alpha * tanh(A*r+W_in*u+bias);
    r_out = r;
    r_out(2:2:end,1) = r_out(2:2:end,1).^2; %even number -> squared
    predict_y = Wout * r_out;
    test_pred_y(t_i,:) = predict_y;
             % u= [predict_y;udata(train_length-1+t_i,41)];
             u= udata(train_length+t_i+1,:);
             u=u';
      % u= predict_y;
%     u=u';
end

error(:,:) = test_pred_y(1:short_prediction_length,:) - test_real_y(1:short_prediction_length,:); 

se_ts = sum( abs(error).^2  ,2);
rmse = sqrt( mean(se_ts) );


testing_start = train_length+valtesting_length+2;

test_pred_yc = zeros(testing_length,dim_T);
test_real_yc = zeros(testing_length,1);
test_real_yc = udata(testing_start:(testing_start+testing_length-1),1:dim_T);

for t_i = 1:testing_length

    r = (1-alpha) * r + alpha * tanh(A*r+W_in*u+bias);
    r_out = r;
    r_out(2:2:end,1) = r_out(2:2:end,1).^2; %even number -> squared
    predict_y = Wout * r_out;
    test_pred_yc(t_i,:) = predict_y;
             u= [predict_y;udata(train_length+valtesting_length-1+t_i,41)];
             
end


%% prediction plots
% 
figure();


target_train=mean(y_train,1);
prediction_train=mean(pred_train,1);

target_test=mean(test_real_y,2);
prediction_test=mean(test_pred_y,2);
% 
target_testc=mean(test_real_yc,2);
prediction_testc=mean(test_pred_yc,2);
c=data(41,:);
plot(c(101:650),target_train,'r')
hold on
plot(c(101:650),prediction_train,'b')
hold on
plot(c(651:700),target_test, 'r--', LineWidth=1);
hold on
plot(c(651:700),prediction_test, 'b--', LineWidth=1);
hold on
plot(c(701:1000),target_testc, 'g-', LineWidth=1);
hold on
plot(c(701:1000),prediction_testc, 'm', LineWidth=1);


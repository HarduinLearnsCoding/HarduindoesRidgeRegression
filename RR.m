%% Solution to HW5 Q1 part d

% By    Anirudh Nakra   UID:118134444 
%       MS,ECE

%% Initialising Variables

clc;
clear;
load('hw05-data.mat');

N=10;
theta=zeros(N,1);
thetaboptimal=zeros(N+1,1);
y=y.';
x=x.';
bigmatrix=zeros(N+1,N+1);
transposed1=ones(N,1);

z=input('Plot for your choice of Gamma,Lambda (Enter 1) or Default choices given in Q1 d (Enter 2)? \n');
resolution=input("Enter the resolution of [-2,2] \n");

switch z
    case 1 %Choice of Gamma and Lambda
        
    gamma=input('Choose kernel parameter gamma\n');
    lambda=input('Choose regularization constant lambda \n');

%% Solving for the optimal theta and bias

[kernelmatrix]=kernelhw05(x,y,N,gamma);
ymatrix=[y;transposed1.'*y];
bigmatrix(1:N,1:N)=kernelmatrix+lambda*eye(N,N);
bigmatrix(1:N,N+1)=transposed1;
bigmatrix(N+1,1:N)=transposed1.'*kernelmatrix;
bigmatrix(N+1,N+1)=N;
thetaboptimal=(bigmatrix)\ymatrix;
b = thetaboptimal(N+1,1);
theta = thetaboptimal(1:N,1);

%% Solving the equations for optimal f(x) and subsequently f(x)+b

xrange=linspace(-2,2,resolution);
kernel=[];

for i=1:N
    for j=1:length(xrange)   
        kernel(j,i)=exp(-gamma*(xrange(1,j)-x(i,1))^2);
    end
end

fx=zeros(resolution,1);

for j=1:length(xrange)
    for i=1:N
        fx(j,1)=theta(i,1)*kernel(j,i)+fx(j,1);
    end
end

for i=1:length(xrange)
     fx(i,1)=fx(i,1)+b;
 end

%% Plotting the relevant variables

figure('units','normalized','outerposition',[0 0 1 1]);

title('Plotting the ridge regression problem');            
h1=plot(x,y,'rh','MarkerSize',15,'MarkerFaceColor','r','DisplayName','Data');
hold on;
h2=plot(xrange,fx,'--bo','MarkerSize',5,'MarkerFaceColor','b','DisplayName','Ridge Regressor');
hold on;
yline(0);
hold off;
xlabel('Data and Optimal Ridge Regressor');
ylabel('Y');
title(['Gamma = ' num2str(gamma)  ' Lambda = ' num2str(lambda)]);
lgndnew=legend([h1(1), h2(1)], 'Data', 'Ridge Regressor');
lgndnew.Position(1) = .8;


case 2    %Prechosen Lambda and Gammas
    
    gamma=[0.1,10,0.1,10];
    lambda=[0.01,1,1,0.01];

%% Solving for the optimal theta and bias

for i=1:length(gamma)

    kernelmatrix(:,:,i)=kernelhw05(x,y,N,gamma(i));
    ymatrix(:,:,i)=[y;transposed1.'*y];
    bigmatrix(1:N,1:N,i)=kernelmatrix(:,:,i)+lambda(i)*eye(N,N);
    bigmatrix(1:N,N+1,i)=transposed1;
    bigmatrix(N+1,1:N,i)=transposed1.'*kernelmatrix(:,:,i);
    bigmatrix(N+1,N+1,i)=N;
    thetaboptimal(:,i)=(bigmatrix(:,:,i))\ymatrix(:,:,i);
    b(i,1) = thetaboptimal(N+1,i);
    theta(:,i) = thetaboptimal(1:N,i);

end
%% Solving the equations for optimal f(x) and subsequently f(x)+b
 
xrange=linspace(-2,2,resolution);
kernel=[];

for z=1:length(gamma)
    for i=1:N
        for j=1:length(xrange)   
            kernel(j,i,z)=exp(-gamma(z)*(xrange(1,j)-x(i,1))^2);
        end
    end
end

fx=zeros(length(xrange),length(gamma));

for z=1:length(gamma)
    for j=1:length(xrange)
        for i=1:N
            fx(j,z)=theta(i,z)*kernel(j,i,z)+fx(j,z);
        end
    end
end

for z=1:length(gamma)
    for i=1:length(xrange)
       fx(i,z)=fx(i,z)+b(z);
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
h1=plot(x,y,'rh','MarkerSize',15,'MarkerFaceColor','r','DisplayName','Data');
hold on;
h2=plot(xrange,fx(:,1),'--bo','MarkerSize',5,'MarkerFaceColor','b','DisplayName','Ridge Regressor');
hold on;
yline(0);
hold off;
xlabel('Data and Optimal Ridge Regressor');
ylabel('Y');
title(['Gamma = ' num2str(gamma(1))  ' Lambda = ' num2str(lambda(1))]);
lgndnew=legend([h1(1), h2(1)], 'Data', 'Ridge Regressor');
lgndnew.Position(1) = .8;


subplot(2,2,2);
h1=plot(x,y,'rh','MarkerSize',15,'MarkerFaceColor','r','DisplayName','Data');
hold on;
h2=plot(xrange,fx(:,2),'--bo','MarkerSize',5,'MarkerFaceColor','b','DisplayName','Ridge Regressor');
hold on;
yline(0);
hold off;
xlabel('Data and Optimal Ridge Regressor');
ylabel('Y');
title(['Gamma = ' num2str(gamma(2))  ' Lambda = ' num2str(lambda(2))]);
lgndnew=legend([h1(1), h2(1)], 'Data', 'Ridge Regressor');
lgndnew.Position(1) = .8;

subplot(2,2,3);
h1=plot(x,y,'rh','MarkerSize',15,'MarkerFaceColor','r','DisplayName','Data');
hold on;
h2=plot(xrange,fx(:,3),'--bo','MarkerSize',5,'MarkerFaceColor','b','DisplayName','Ridge Regressor');
hold on;
yline(0);
hold off;
xlabel('Data and Optimal Ridge Regressor');
ylabel('Y');
title(['Gamma = ' num2str(gamma(3))  ' Lambda = ' num2str(lambda(3))]);
lgndnew=legend([h1(1), h2(1)], 'Data', 'Ridge Regressor');
lgndnew.Position(1) = .8;

subplot(2,2,4);
h1=plot(x,y,'rh','MarkerSize',15,'MarkerFaceColor','r','DisplayName','Data');
hold on;
h2=plot(xrange,fx(:,4),'--bo','MarkerSize',5,'MarkerFaceColor','b','DisplayName','Ridge Regressor');
hold on;
yline(0);
hold off;
xlabel('Data and Optimal Ridge Regressor');
ylabel('Y');
title(['Gamma = ' num2str(gamma(4))  ' Lambda = ' num2str(lambda(4))]);
lgndnew=legend([h1(1), h2(1)], 'Data', 'Ridge Regressor');
lgndnew.Position(1) = .8;
end






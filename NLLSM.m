function a = NLLSM(database) 
% An iterative method to solve non-linear LSM problem

% a: the last value vector for [m, l, beta] when iteration stops

% RUN: NLLSM("LSM_data.mat")

clear all
close all
% load db with fields x and y
db = load("LSM_data.mat");

% TODO: x: matrix, each column vector contains [q, qp, qpp] 
% from an experiment with pendulum dynamics model
x = db.x;
% TODO: y: column vector as input, it contains real tau values 
% from an experiment with pendulum dynamics model
y = db.y;

tol = 0.00001;  %--set a tolerance value for the accuracy, DO NOT CHANGE IT

% TODO: Try with different learning rate lamda, e.g., 0.5, 0.2, 0.1, 0.05, 0.001 ...
lamda = 0.1;  %--set a learning rate for incrementing 'a'

% TODO: set initial guess for a=[m, beta, l]. You're free to set a new inital
% guess. The initial values must agree with the
% variable definitions. What is the physical meaning of each variable?
a = [10 0.5 10];   


% TODO: Try with different iteration numbers iter_max, e.g., 30, 40, 50 ,70, 80, 100, 150 so on ...
iter_max = 175;  %--set maximum iteration number to run for
n = length(x);  %--number of data samples, each sample contains [q, qp, qpp]   

% Fixed parameters, do not change!
I = 0.04;  %-- inertial coefficient
g = 9.81;  %-- gravitational acceleration

% sampling time 1Khz (in seconds)
Dt=0.001;


% TODO: Extract the joint positions, velocitites, and accelerations from
% the matrix x=[q qp qpp]\in (nx3)
q = x(:,1);  %-- joint position
qp = x(:,2);  %-- joint velocity
qpp = x(:,3);  %-- joint acceleration

% save each parameter iteration and the estimation error for plotting
va(1,:)=a;
ve=[];
%% starting the iteration process
for iter = 1:iter_max
    
    % TODO: Assign the elements of the vector a to the robot parametes m,beta,l.
    m = a(1);  %-- mass 
    beta = a(2);  %-- viscous friction coefficient 
    l = a(3);  %-- length 
    
    for i = 1:n
	% TODO: predicted tau value (\hat y), you need 'qpp(i)', 'qp(i)' and 'q(i)'.
	% Use the model from Eq. 2.1
        f(i) = (m*l^2+I)*qpp(i)+beta*qp(i)+m*g*l*sin(q(i));
        
	%% A should be a n x 3 matrix, n stands for the number of samples, 
        % each row contains a partial derivative of tau w.r.t --> m, l and beta
       
        % TODO: calculate partial derivative of tau w.r.t. 'm'. The partial
        % derivative will be a function of qpp(i) and q(i)
		A(i, 1) = l^2*qpp(i)+g*l*sin(q(i));
		
		% TODO: calculate partial derivative of tau w.r.t. 'beta'. The partial
        % derivative will be a function of qp(i) 
        A(i, 2) = qp(i);
		
		% TODO: calculate partial derivative of tau w.r.t. 'l'. The partial
        % derivative will be a function of qpp(i) and q(i)
        A(i, 3) = 2*m*l*qpp(i)+m*g*sin(q(i));
		
		% TODO: calculate the error (residual) given real tau value 'y(i)' and 
        % the predicted tau value '\hat y=f(i)'
        e(i) = y(i)-f(i);
    end
	
    % Save the vector of errors to plot them
    ve=[ve e];

    % TODO: calculate the increment 'Delta_a' for 'a': [m, l, beta] using
    % the matrix 'A' and the error 'e' with the LSM
    Delta_a = inv(A'*A)*A'*e'; 
    % display 
    fprintf(1,'iteration: %d\n', iter);
    fprintf(1,'m, l, beta: \n');
    disp(a)
    fprintf(1,'increment for a: \n');
    disp(Delta_a')
	
    % TODO: calculate the new approximation for 'a', using 'Delta_a' and
    % the learning factor lambda
    a = a +lamda*Delta_a';
    
    % Save the last vector 'a' for plotting
    va(iter+1,:)=a;
	
    % If the increments of 'a' ('Delta_a') are less than a tolerance break
    % the iteration
    % If the algorithm reaches iter_max without reaching the tolerance
    % value, then, your algorithm is not converging well
 
    if (abs(Delta_a(1)) < tol && abs(Delta_a(2)) < tol && abs(Delta_a(3)) < tol)
        disp("Iteration has converged")
        break
        
    end
end

%% Plotting results (Nothing to do here!)
% NOTE: if you get an error here, check the dimension of  y, q, qp, and qpp
% Also, if the plots look out of scale, check the value of DT, 
% it must be 1KHz in seconds.

qe_old=q(1);
qep_old=qp(1);
qepp_old=qpp(1);

for i=1:n
    
    % Time vector
    t(i)=(i-1)*Dt;
    
    % 
    qepp=(1/(I+m*l^2))*(y(i)-beta*qep_old-m*l*g*sin(qe_old));
    
    qep=qep_old+(qepp+qepp_old)*Dt/2.0;
    
    qe=qe_old+(qep+qep_old)*Dt/2.0;
    
    Qpp(i)=qepp;
    Qp(i)=qep;
    Q(i)=qe;
    
    qe_old=qe;
    qep_old=qep;
    qepp_old=qepp;
    
    eq(i)=q(i)-qe;
    eqp(i)=qp(i)-qep;
    eqpp(i)=qpp(i)-qepp;

end

tfontSize=20;
lfontSize=16;
lwidth=2;

figure(3)

plot(ve, 'b -','LineWidth',lwidth)
grid on
title('Estimation error convergency','FontSize', tfontSize);
xlabel('time [s]','FontSize', lfontSize)
ylabel('e=y-f [Nm]','Interpreter','latex','FontSize', lfontSize)
legend('e','FontSize', lfontSize);


figure(1)

subplot(4,1,1); plot(t,Q,'b -',t,x(:,1),'r --','LineWidth',lwidth)
grid on
title('Joint Position $\hat{q}$ vs q','Interpreter','latex','FontSize', tfontSize);
xlabel('time [s]','FontSize', lfontSize)
ylabel('$q$ [rad]','Interpreter','latex','FontSize', lfontSize)
legend('$\hat{q}$','q','Interpreter','latex','FontSize', lfontSize);

subplot(4,1,2); plot(t,Qp,'b -',t,x(:,2),'r --','LineWidth',lwidth)
grid on
title('Joint Velocity $\hat{\dot{q}}$ vs $\dot{q}$','Interpreter','latex','FontSize', tfontSize);
xlabel('time [s]','FontSize', lfontSize)
ylabel('$\hat{\dot{q}}$ [rad/s]','Interpreter','latex','FontSize', lfontSize)
legend('$\hat{\dot{q}}$','$\dot q$','Interpreter','latex','FontSize', lfontSize);


subplot(4,1,3); plot(t,Qpp,'b -',t,x(:,3),'r --','LineWidth',lwidth)
grid on
title('Joint Acceleration $\hat{\ddot{q}}$ vs $\ddot{q}$','Interpreter','latex','FontSize', tfontSize);
xlabel('time [s]','FontSize', lfontSize)
ylabel('$\hat{\ddot{q}}$ [rad/s]','Interpreter','latex','FontSize', lfontSize)
legend('$\hat{\ddot{q}}$','$\ddot q$','Interpreter','latex','FontSize', lfontSize);

subplot(4,1,4);plot(t,rad2deg(eq),'b -',t,rad2deg(eqp),'r -', t,rad2deg(eqpp),'k -','LineWidth',lwidth)
title('Estimation Errors','Interpreter','latex','FontSize', tfontSize)
xlabel('time [s]','FontSize', lfontSize);
ylabel('error $[deg, \frac{deg}{s}, \frac{deg}{s^2}]$','Interpreter','latex','FontSize', lfontSize)
legend('$e$','$\dot e$','$\ddot e$','Interpreter','latex','FontSize', lfontSize);
grid on


%% Plotting parameter convergence
figure(2)

plot3(va(1,1),va(1,2),va(1,3), 'Marker','o',...
            'Color',[0 1 0.3],'MarkerSize',10,'LineWidth',4);
hold on
plot3(va(2:end-1,1),va(2:end-1,2),va(2:end-1,3), 'Marker','+',...
            'Color',[1,0.3,0],'MarkerSize',10,'LineWidth',2);
plot3(va(end,1),va(end,2),va(end,3), 'Marker','d',...
            'Color',[0 0 0],'MarkerSize',10,'LineWidth',4);
grid on
xlabel('m');
ylabel('$\beta$','Interpreter','latex');
zlabel('l');


end
    

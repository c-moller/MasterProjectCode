tic

clear
tps =  20;


P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',25, 'k_cs', 5, 'k_cf', 25,...
    'RBS',5,'L', 50, 'epsi', 10.^-2);

%Initial values 
y0 = [0 2 10 0];

%ODE solver
[T1, x] = ode23(@(t,y) ODEdistRibo(t, y, P1), [0 tps], y0);

%Parameter Extract

for i = 1: length(T1)
[x2, paraout] = ODEdistRibo(T1(i),x(i,:),P1);
para_init(i,1) = paraout(1);
para_release(i,1) = paraout(2);
para_RBSr(i,1) = paraout(3);
end
 
%Plots
figure
plot (T1, para_init, 'k'), hold on
plot (T1, para_release, 'g')
plot (T1, x(:, 3), 'r')

xlabel ('Time');
ylabel ('[]');
legend ('Initiation Rate','Release Rate', 'TC');


clear 
P2 = struct('tps', 20, 'var', 20);
Transient_Sweep(P2);
RateRBS_Sweep(P2);



toc
tic
%close all;
clear;

tps =  150;

%Parameters

P1 = struct('ngenes', 4,'s0', 5, 'mRNAPool', [1 1],...
	'x_cs',[50 50], 'k_cs', [1 1], 'k_cf', [1 1],...
	'RBS', [2 10], 'L', [100 100],...  
    'vs', 2, 'Ks', 1, 'tag', 1, 'Kdeg', 10,'CellMass', 100);

%Initial values 
%y0 = ones(1,(1+2*P1.ngenes));
y0 = [100 5 1 1 1 1];

[T1, x] = ode15s(@(t,y) ODE1(t, y, P1), [0 tps],y0);



for i = 1: length(T1)
[x2, paraout] = ODE1(T1(i),x(i,:),P1);
para_init(i,1) = paraout(1,1);
para_init(i,2) = paraout(1,2);
para_tpsmRNA(i,1) = paraout(2,1);
para_tpsmRNA(i,2) = paraout(2,2);
para_elon(i,1) = paraout(3,1);
para_elon(i,2) = paraout(3,2);
para_Queue(i,1) = paraout(4,1);
para_Queue(i,2) = paraout(4,2);
end

figure 
subplot(2,2,1)
plot(T1,para_init);
legend ('Meta','Syn');
xlabel ('Time');
ylabel ('Initiation rate');
subplot(2,2,2)
plot(T1,para_tpsmRNA);
legend ('Meta','Syn');
xlabel ('Time');
ylabel ('mRNA translation time');
subplot(2,2,3)
plot(T1,para_elon);
legend ('Meta','Syn');
xlabel ('Time');
ylabel ('Elongation rate');
subplot(2,2,4)
plot(T1,para_Queue);
legend ('Meta','Syn');
xlabel ('Time');
ylabel ('Ribosome on mRNA');

%%
%Analysis

%Plots
figure
subplot(2,2,1)
plot (T1, x(:, 2))
legend ('Ribosomes');
xlabel ('Time');
ylabel ('[]');
subplot(2,2,2)
plot (T1, x(:, 3)), hold on
plot (T1, x(:, 5))
%plot (T1, x(:, 6))
%title ('Translation complexes and ribosome concentration');
xlabel ('Time');
ylabel ('[TC]');
legend ( 'Meta TC', 'Synth. TC'); %, 'Ribo TC'

subplot(2,2,3)
plot (T1, x(:, 1))
xlabel ('Time');
ylabel ('[]');
legend ('Amino Acid');
subplot(2,2,4)
plot (T1, x(:, 4)), hold on
plot (T1, x(:, 6))
%title ('Amino Acids, Metabolic, Synth. proteins concentrations (MS)');
xlabel ('Time');
ylabel ('[]');
legend ('Meta', 'synth');

toc
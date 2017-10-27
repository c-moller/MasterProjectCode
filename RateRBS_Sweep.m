function RateRBS_Sweep(P2)

tps=1000;
y = 0:0.1:10;
endSweep = length(y);
j=1;

for var = 1:endSweep
    
P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',25, 'k_cs', 5, 'k_cf', 25,...
    'RBS',y(j),'L', 50, 'epsi', 1e-2, 'km',5);
    
    %Initial values
    y0 = [P1.L 1 0 0];
    
    %ODE solver
    %[T1, x] = ode23(@(t,y) ODEdistRibo(t, y, P1), [0 tps], y0);
    [T1, x] = ode23(@(t,y) ProteaseODE(t, y, P1), [0 tps], y0);
    
    %Parameter Extract
    
    for i = 1: length(T1)
        %[x2, paraout] = ODEdistRibo(T1(i),x(i,:),P1);
        [x2, paraout] = ProteaseODE(T1(i),x(i,:),P1);
        Time(i,j)=T1(i);
        para_init(i,j) = paraout(1);
        para_release(i,j) = paraout(2);
        para_RBSr(i,j) = paraout(3);
    end
    last_TC(j,:)= x(i,3);
    last_init(j,:)= para_init(i,j);
    last_release(j,:)= para_release(i,j);
    last_RBSr(j,:)= para_RBSr(i,j);
    
    j=j+1;
end
figure
    xplot =last_RBSr;
    y1plot = last_TC;
    y2plot = last_init;
    y3plot =  last_release;
    subplot(3,1,1);
    plot(xplot,y1plot, '+ r'),hold on
    xlabel('RBS Binding Rate (k_{RBS}*Ribo_{free})');
    ylabel('TC');
    title({'k_{cs} =5 , k_{cf} =25, mRNA length =50, Slow Codon Location =25';'Initial values:AA Pool=50, Ribo Pool=1, TC=0'});
    subplot(3,1,2);
    plot(xplot,y2plot, '+ b'),hold on
    xlabel('RBS Binding Rate (k_{RBS}*Ribo_{free})');
    ylabel('Initiation Rate');
    subplot(3,1,3);
    plot(xplot,y3plot, '+ g '),hold on
    xlabel('RBS Binding Rate (k_{RBS}*Ribo_{free})');
    ylabel('Release Rate');



end





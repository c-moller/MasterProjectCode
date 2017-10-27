function [xdot, outpara] = ODEstepMin(~, y, P1)

Amino = y(1);              %Amino acids
Ribo = y(2);  %Ribosomes
TC = y(3); %TC Meta.
Meta = y(4);  %Metabolic Proteins
mRNA = P1.mRNAPool;  %mRNA Metabolic Proteins

%Rates
ribo_density = TC./P1.mRNAPool;
ribo_density_max = P1.x_cs; %Max queue length, assume length of mRNA in codon increments so L = 100 is mRNA can have 100TC on it
% queue_function = ribo_density_max^P1.n ./(ribo_density_max^P1.n + ribo_density^P1.n);    %is 1 when no queue, n steepness => higher value queue close to 0
% k_init = P1.RBS .* Ribo .* P1.mRNAPool .* queue_function; %Ribosome initiation rate
queue_function=1;

% RBS = (ribo_density < P1.x_cs)*P1.RBS
k_init = (ribo_density < P1.x_cs)*P1.RBS*Ribo;
k_trans_slow =P1.k_cs;% * Amino;
tps_cs = 1/(P1.k_cs);%Time for single ribosome to translate slow codon
tps_cf = 1/(P1.k_cf); %Time for single ribosome to translate fast codon
%k_init = min(k_RBS, (k_trans_slow/mRNA)); % Initiation per mRNA
k_elon = min(k_init, (k_trans_slow/mRNA));
k_Prod_Meta = k_elon*TC;

% queue_function = (ribo_density < P1.x_cs)*P1.x_cs + 1;
% tps_avgmRNA = (P1.L-queue_function)*tps_cf+queue_function*tps_cs;
% k_elon = 1/tps_avgmRNA;
% k_Prod_Meta = k_elon*TC;

k_RBS=P1.RBS*Ribo;
% queue_function = 0;
% ribo_density_max = 6;


%-------------------------------------------------------------------------------------------------------------------------------------------------
dAmino =  0; %  k_Acat - dilution*Amino;
dRibo = 0;%k_Prod_Meta - k_init*mRNA_Meta - dilution*Ribo; %
dTC_Meta =  k_init - k_Prod_Meta;% - dilution*TC_Meta;
dMeta =  k_Prod_Meta;% - dilution*Meta;


%-------------------------------------------------------------------------------------------------------------------------------------------------

xdot = [dAmino; %x(1)
    dRibo;%x(2)
    dTC_Meta;%x(3)
    dMeta; %x(4)
    ];
outpara = [k_init;
    k_elon;
    queue_function;
    ribo_density;
    k_RBS;
    ribo_density_max;
    ];

end


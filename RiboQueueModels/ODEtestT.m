function [xdot, outpara] = ODEtestT(~, y, P1)

Amino = y(1);              %Amino acids
Ribo = y(2);  %Ribosomes
TC = y(3); %TC Meta.
Meta = y(4);  %Metabolic Proteins
mRNA_Meta = P1.mRNAPool;  %mRNA Metabolic Proteins



%Rates
ribo_density = TC./P1.mRNAPool;
RBS_Meta = (ribo_density < P1.x_cs)*P1.RBS;
tps_cs = 1/(P1.k_cs);%Time for single ribosome to translate slow codon
tps_cf = 1/(P1.k_cf); %Time for single ribosome to translate fast codon
k_init = RBS_Meta*Ribo; % Initiation per mRNA

queue_function = (k_init > P1.k_cs)*P1.x_cs + 1;
tps_avgmRNA = (P1.L-queue_function)*tps_cf+queue_function*tps_cs;
k_elon = 1/tps_avgmRNA;
k_Prod_Meta = k_elon*TC;

RiboRBS = P1.RBS*Ribo;
ribo_density_max = (1/tps_cs);

%-------------------------------------------------------------------------------------------------------------------------------------------------
dAmino =  0; %  k_Acat - dilution*Amino;
dRibo = 0;%k_Prod_Meta - k_init*mRNA_Meta - dilution*Ribo; %
dTC_Meta =  k_init*mRNA_Meta - k_Prod_Meta;% - dilution*TC_Meta;
dMeta =  k_Prod_Meta;% - dilution*Meta;

%-------------------------------------------------------------------------------------------------------------------------------------------------

xdot = [dAmino; %x(1)
    dRibo;%x(2)
    dTC_Meta;%x(3)
    dMeta; %x(4)
    ];
outpara = [k_init;
    k_Prod_Meta;
    queue_function;
    ribo_density;
    RiboRBS;
    ribo_density_max;
    ];

end


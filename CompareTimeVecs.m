%compare results from time vectors with 100, 50, and 25 points for all five
%distributions in the zero noise case 
%this is with the original forward problem DE



%Create original distribution


%Bigaussian
gls_optpar_bi_100 =[0.000000000000000,
   0.006178125635662
   0.150004317399579
   0.122031536345074
   0.088607951556955
   0.097569446496761
   0.129164072279879
   0.151319297782558
   0.145558848488781
   0.105200019855914
   0.004265092517468
   0.000101291641369]

   gls_optpar_bi_50 = [      0
   0.000000000000000
   0.098188881191274
   0.106286933433693
   0.091124627637342
   0.085442786230165
   0.094407545587683
   0.109294397120783
   0.119615975570285
   0.118889250641031
   0.104610122291207
   0.072139480296536
   0.000000000000000
   0.000000000000000]

    gls_optpar_bi_25=[0
   0.112157388072792
   0.178970665389985
   0.204596422946065
   0.205985680132871
   0.194784696338348
   0.102860296240417
   0.000644850879521 ]        

[sprobs_bi_t100] = DistFn2('Bigaussian',linspace(0,1,101),0,1); %this now has nothing to do with t100 because I changed the number of points. 
[sprobs_n_t100] = DistFn2('Normal',linspace(0,1,101),0,1);
[sprobs_u_t100] = DistFn2('Uniform',linspace(0,1,101),0,1);

%Normal
gls_optpar_n_100 = [0.000315778838166
                   0
                   0
                   0
                   0
   0.033368659192725
   0.133026335571202
   0.207167894893708
   0.239721558744881
   0.220653713707274
   0.145836528877979
   0.016175065913135
   0.003734238194837
   0.000000000000000
   0.000000000000000
                   0
   0.000000226066093]

   gls_optpar_n_50 = [0.000325388582934
                   0
                   0
                   0
   0.083719939922603
   0.249091489374704
   0.323219623875663
   0.268853541914245
   0.074790016329850
                   0
                   0
                   0
   0.000000000000000]

   gls_optpar_n_25 =  [0.000327290102750
   0.000000090337830
   0.000000046263959
   0.000000567012114
   0.083664799801780
   0.249186660946973
   0.323212225565829
   0.268785211009375
   0.074821363123871
   0.000001745835521
   0.000000000000000
                   0
   0.000000000000000]

%OnePoint
gls_optpar_one_100 =  [0
                   0
                   0
                   0
                   0
                   0
   0.000000000000000
                   0
   0.029297244476515
   0.970702755523485
                   0
   0.000000000000000
                   0
                   0
                   0
   0.000000000000000
                   0
   0.000000000000000
                   0
                   0
                   0
                   0
   0.000000000000000
                   0]

   gls_optpar_one_50 = [0
   0.000000000000000
                   0
   0.000000000000000
   0.000000000000000
   0.000000000000000
                   0
                   0
   0.029297657928656
   0.970702342071344
   0.000000000000000
   0.000000000000000
                   0
                   0
                   0
                   0
                   0
                   0
   0.000000000000000
   0.000000000000000
   0.000000000000000
   0.000000000000000
   0.000000000000000
                   0]                

                   gls_optpar_one_25 =   [0
                   0
   0.000000000000000
   0.000000000000000
   0.000000000000000
                   0
                   0
                   0
   0.029299326133966
   0.970700673866034
                   0
                   0
   0.000000000000000
   0.000000000000000
   0.000000000000000
                   0
   0.000000000000000
   0.000000000000000
                   0
   0.000000000000000
   0.000000000000000
   0.000000000000000
                   0
   0.000000000000000]
%TwoPoints

gls_optpar_two_100 =   [0.000000000000000
   0.005735359774547
   0.338204224192827
   0.054868767485045
   0.052093546250630
   0.218097095396852
   0.331001006900099
   0.000000000000000]

   gls_optpar_two_50 =  [0
   0.006570680549118
   0.334578730604113
   0.059330685728350
   0.052598050300170
   0.213972652461578
   0.332949200356670
                   0]
gls_optpar_two_25 = [0
   0.006418950601825
   0.335801931884689
   0.057193688669505
   0.053159930506698
   0.215113736075157
   0.332311762262126
                   0]
%Uniform
gls_optpar_u_100 = [ 0.025216360354617
   0.033059228778641
   0.037709986379885
   0.039948410141271
   0.040511947295979
   0.040053460287159
   0.039109555591315
   0.038082846389527
   0.037239356223314
   0.036720001099160
   0.036562767544498
   0.036731192773684
   0.037144050228966
   0.037702005120481
   0.038308278737535
   0.038881985812021
   0.039364211058362
   0.039718189395403
   0.039925565307946
   0.039662635554869
   0.038832880981063
   0.037913190120325
   0.036894913378979
   0.035762293717335
   0.034491488028655
   0.033050871427817
   0.031402328271193]

   gls_optpar_u_50 =  [0.035681255596575
   0.077701623747903
   0.085583568419327
   0.079219587741045
   0.072085097503793
   0.070000694887597
   0.072882319334393
   0.077932867839650
   0.082351181288823
   0.082569039771989
   0.078251276080952
   0.071806936617920
   0.062949628840074
   0.050984922329957]

   gls_optpar_u_25 =  [0.038929327498944
   0.098083718324332
   0.100503795707939
   0.087476665882440
   0.081317061702648
   0.085696999315951
   0.094499564468312
   0.097722131696090
   0.095785144662026
   0.088383606383979
   0.075474909783951
   0.056127074573390]

   % %bilist = (gls_optpar_bi_100, gls_optpar_bi_50, gls_optpar_bi_25)
   % for gls_optpar_bi_100, gls_optpar_bi_50, gls_optpar_bi_25
   %     1+1
   % end

   figure
   yyaxis left
   stem(linspace(0,1,length(gls_optpar_bi_100)),gls_optpar_bi_100*length(gls_optpar_bi_100),'--db','MarkerSize',10,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_bi_50)),gls_optpar_bi_50*length(gls_optpar_bi_50),'--*b','MarkerSize',8,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_bi_25)),gls_optpar_bi_25*length(gls_optpar_bi_25),'--ob','MarkerSize',7,'LineWidth',2)
     hold on
     yyaxis right
     plot(linspace(0,1,101),sprobs_bi_t100,'r')
     legend('100 time points','50 time points','25 time points')
 

     figure
     yyaxis left
   stem(linspace(0,1,length(gls_optpar_n_100)),gls_optpar_n_100*length(gls_optpar_n_100),'--d','MarkerSize',10,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_n_50)),gls_optpar_n_50*length(gls_optpar_n_50),'--*','MarkerSize',8,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_n_25)),gls_optpar_n_25*length(gls_optpar_n_25),'--o','MarkerSize',7,'LineWidth',2)
     hold on
     ylabel('Recovered Prop. of Pop.')
     yyaxis right
     plot(linspace(0,1,101),sprobs_n_t100,'k')
     legend('100 time points','50 time points','25 time points')
     title('Recovered Normal PMF for Different Time Steps')
     xlabel('Sensitivity to Treatment {\it s}')
     ylabel('Original Prop. of Pop.')

     figure
   stem(linspace(0,1,length(gls_optpar_one_100)),gls_optpar_one_100,'--d','MarkerSize',10,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_one_50)),gls_optpar_one_50,'--*','MarkerSize',8,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_one_25)),gls_optpar_one_25,'--o','MarkerSize',7,'LineWidth',2)
     legend('100 time points','50 time points','25 time points')

       figure
   stem(linspace(0,1,length(gls_optpar_two_100)),gls_optpar_two_100,'--d','MarkerSize',10,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_two_50)),gls_optpar_two_50,'--*','MarkerSize',8,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_two_25)),gls_optpar_two_25,'--o','MarkerSize',7,'LineWidth',2)
     legend('100 time points','50 time points','25 time points')

            figure
            yyaxis left
   stem(linspace(0,1,length(gls_optpar_u_100)),gls_optpar_u_100*length(gls_optpar_u_100),'--d','MarkerSize',10,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_u_50)),gls_optpar_u_50*length(gls_optpar_u_50),'--*','MarkerSize',8,'LineWidth',2)
   hold on
     stem(linspace(0,1,length(gls_optpar_u_25)),gls_optpar_u_25*length(gls_optpar_u_25),'--o','MarkerSize',7,'LineWidth',2)
     ylim([0 4*median(gls_optpar_u_100)*length(gls_optpar_u_100)])  
     yyaxis right
     plot(linspace(0,1,101),sprobs_u_t100,'k')
     ylim([0 4*median(sprobs_u_t100)])
     legend('100 time points','50 time points','25 time points')
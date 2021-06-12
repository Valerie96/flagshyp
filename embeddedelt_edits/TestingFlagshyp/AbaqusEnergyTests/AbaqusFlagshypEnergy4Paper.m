%%  Strain Rate
set(0,'defaultfigurecolor',[1 1 1]);

Abq1h_R1 = ReadHost2("Ab_1h_Rate1");
Abq1h_R5 = ReadHost2("Ab_1h_OG");
Abq1h_R10 = ReadHost2("Ab_1h_Rate10");
Abq1h_R25 = ReadHost2("Ab_1h_Rate25");
Abq1h_R100 = ReadHost2("Ab_1h_Rate100");
Abq1h_R200 = ReadHost2("Ab_1h_Rate200");
graphsize=[100 100 800 400];
name1 = "1h Rate1";
name2 = "1h Rate5";
name3 = "1h Rate10";
name4 = "1h Rate25";

% PlotEnergy3([Abq1h_R1.time, Abq1h_R1.KE], [Abq1h_R5.time, Abq1h_R5.KE], [Abq1h_R10.time, Abq1h_R10.KE], name1, name2,name3,'Kinetic Energy')
% PlotEnergy3([Abq1h_R1.time, Abq1h_R1.IE], [Abq1h_R5.time, Abq1h_R5.IE], [Abq1h_R10.time, Abq1h_R10.IE], name1, name2,name3,'Internal Energy')
% PlotEnergy3([Abq1h_R1.time, Abq1h_R1.WK], [Abq1h_R5.time, Abq1h_R5.WK], [Abq1h_R10.time, Abq1h_R10.WK], name1, name2,name3,'External Energy')
% PlotEnergy3([Abq1h_R1.time, Abq1h_R1.ETOTAL], [Abq1h_R5.time, Abq1h_R5.ETOTAL], [Abq1h_R10.time, Abq1h_R10.ETOTAL], name1, name2,name3,'Total Energy')

[~, ~, Abq1h_25t_R1] = ReadHostTruss2("Ab_1h_25t_Rate1");
[~, ~, Abq1h_25t_R5] = ReadHostTruss2("Ab_1h_25t");
[~, ~, Abq1h_25t_R10] = ReadHostTruss2("Ab_1h_25t_Rate10");
[~, ~, Abq1h_25t_R25] = ReadHostTruss2("Ab_1h_25t_Rate25");
[~, ~, Abq1h_25t_R100] = ReadHostTruss2("Ab_1h_25t_Rate100");
[~, ~, Abq1h_25t_R200] = ReadHostTruss2("Ab_1h_25t_Rate200");
graphsize=[100 100 800 400];
name5 = "1h 25t Rate1";
name6 = "1h 25t Rate5";
name7 = "1h 25t Rate10";
name8 = "1h 25t Rate25";

PlotEnergy3([Abq1h_25t_R10.time, Abq1h_25t_R10.KE], [Abq1h_25t_R100.time, Abq1h_25t_R100.KE], [Abq1h_25t_R200.time, Abq1h_25t_R200.KE], name5, name6,name7,'Kinetic Energy')
% PlotEnergy3([Abq1h_25t_R1.time, Abq1h_25t_R1.IE], [Abq1h_25t_R5.time, Abq1h_25t_R5.IE], [Abq1h_25t_R10.time, Abq1h_25t_R10.IE], name5, name6,name7,'Internal Energy')
% PlotEnergy3([Abq1h_25t_R1.time, Abq1h_25t_R1.WK], [Abq1h_25t_R5.time, Abq1h_25t_R5.WK], [Abq1h_25t_R10.time, Abq1h_25t_R10.WK], name5, name6,name7,'External Energy')
% PlotEnergy3([Abq1h_25t_R1.time, Abq1h_25t_R1.ETOTAL], [Abq1h_25t_R5.time, Abq1h_25t_R5.ETOTAL], [Abq1h_25t_R10.time, Abq1h_25t_R10.ETOTAL], name5, name6,name7,'Total Energy')

%%
    AveKE_0t = Abq1h_R1.KE(end);  AveKE1_25t = Abq1h_25t_R1.KE(end);
    AveKE5_0t = Abq1h_R5.KE(end);  AveKE5_25t = Abq1h_25t_R5.KE(end);
    AveKE10_0t = Abq1h_R10.KE(end);AveKE10_25t = Abq1h_25t_R10.KE(end);
    AveKE25_0t = Abq1h_R25.KE(end); AveKE25_25t = Abq1h_25t_R25.KE(end);
    AveKE100_0t = Abq1h_R100.KE(end);AveKE100_25t = Abq1h_25t_R100.KE(end);
    AveKE200_0t = Abq1h_R200.KE(end);AveKE200_25t = Abq1h_25t_R200.KE(end);
    
    figure; hold on; grid on;
    plot(1, Abq1h_R1.IE(end),'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(5, Abq1h_R5.IE(end),'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(10, Abq1h_R10.IE(end),'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(25, Abq1h_R25.IE(end),'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(100, Abq1h_R100.IE(end),'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(200, Abq1h_R200.IE(end),'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    
    plot(1, Abq1h_25t_R1.IE(end),'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(5, Abq1h_25t_R5.IE(end),'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(10, Abq1h_25t_R10.IE(end),'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(25, Abq1h_25t_R25.IE(end),'r.','MarkerSize',25, 'DisplayName', '25 Fibers');
    plot(100, Abq1h_25t_R100.IE(end),'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(200, Abq1h_25t_R200.IE(end),'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    xlabel("Strain Rate (1/s)");
    ylabel("Energy (J)");
    title("Internal Energy vs Applied Strain Rate");
    legend('show');
    
    figure; hold on; grid on;
    plot(1, AveKE_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(5, AveKE5_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(10, AveKE10_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(25, AveKE25_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(100, AveKE100_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(200, AveKE200_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    
    plot(1, AveKE1_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(5, AveKE5_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(10, AveKE10_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(25, AveKE25_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers');
    plot(100, AveKE100_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(200, AveKE200_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    xlabel("Strain Rate (1/s)");
    ylabel("Energy (J)");
    title("Kinetic Energy vs Applied Strain Rate");
    legend('show');

    %% Kinetic Energy Rate
    %Energy Difference (Error)
    AveTotEnergy_0t = Abq1h_R1.KE(end);  AveTotEnergy1_25t = Abq1h_25t_R1.KE(end);
    AveTotEnergy5_0t = Abq1h_R5.KE(end);  AveTotEnergy5_25t = Abq1h_25t_R5.KE(end);
    AveTotEnergy10_0t = Abq1h_R10.KE(end);AveTotEnergy10_25t = Abq1h_25t_R10.KE(end);
    AveTotEnergy25_0t = Abq1h_R25.KE(end); AveTotEnergy25_25t = Abq1h_25t_R25.KE(end);
    AveTotEnergy100_0t = Abq1h_R100.KE(end);AveTotEnergy100_25t = Abq1h_25t_R100.KE(end);
    AveTotEnergy200_0t = Abq1h_R200.KE(end);AveTotEnergy200_25t = Abq1h_25t_R200.KE(end);
    
    AveErr1 = abs(AveTotEnergy1_25t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr5 = abs(AveTotEnergy5_25t-AveTotEnergy5_0t)/abs(AveTotEnergy5_0t);
    AveErr10 = abs(AveTotEnergy10_25t-AveTotEnergy10_0t)/abs(AveTotEnergy10_0t);
    AveErr25 = abs(AveTotEnergy25_25t-AveTotEnergy25_0t)/abs(AveTotEnergy25_0t);
    AveErr100 = abs(AveTotEnergy100_25t-AveTotEnergy100_0t)/abs(AveTotEnergy100_0t);
    AveErr200 = abs(AveTotEnergy200_25t-AveTotEnergy200_0t)/abs(AveTotEnergy200_0t);

    figure; hold on; grid on;
    plot(1, AveErr1*100,'.','MarkerSize',25, 'DisplayName', '1 1/s'); 
    plot(5, AveErr5*100,'.','MarkerSize',25, 'DisplayName', '5 1/s'); 
    plot(10, AveErr10*100,'.','MarkerSize',25, 'DisplayName', '10 1/s'); 
    plot(25, AveErr25*100,'.','MarkerSize',25, 'DisplayName', '25 1/s'); 
    plot(100, AveErr100*100,'.','MarkerSize',25, 'DisplayName', '100 1/s'); 
    plot(200, AveErr200*100,'.','MarkerSize',25, 'DisplayName', '200 1/s');     
    xlabel("Strain Rate (1/s)");
    ylabel("Difference (Error) in Energy (%)");
    title("Difference in Kinetic Energy vs Applied Strain Rate");
    legend('show');

    figure; hold on; grid on;
    plot(1, AveTotEnergy_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(5, AveTotEnergy5_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(10, AveTotEnergy10_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(25, AveTotEnergy25_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(100, AveTotEnergy100_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(200, AveTotEnergy200_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    
    plot(1, AveTotEnergy1_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(5, AveTotEnergy5_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(10, AveTotEnergy10_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(25, AveTotEnergy25_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers');
    plot(100, AveTotEnergy100_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(200, AveTotEnergy200_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    xlabel("Strain Rate (1/s)");
    ylabel("Energy (J)");
    title("Kinetic Energy vs Applied Strain Rate");
    legend('show');
    
    
%% Internal Energy Rate
    %Energy Difference (Error)
    AveTotEnergy_0t = Abq1h_R1.IE(end);  AveTotEnergy1_25t = Abq1h_25t_R1.IE(end);
    AveTotEnergy5_0t = Abq1h_R5.IE(end);  AveTotEnergy5_25t = Abq1h_25t_R5.IE(end);
    AveTotEnergy10_0t = Abq1h_R10.IE(end);AveTotEnergy10_25t = Abq1h_25t_R10.IE(end);
    AveTotEnergy25_0t = Abq1h_R25.IE(end); AveTotEnergy25_25t = Abq1h_25t_R25.IE(end);
    AveTotEnergy100_0t = Abq1h_R100.IE(end);AveTotEnergy100_25t = Abq1h_25t_R100.IE(end);
    AveTotEnergy200_0t = Abq1h_R200.IE(end);AveTotEnergy200_25t = Abq1h_25t_R200.IE(end);
    
    AveErr1 = abs(AveTotEnergy1_25t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr5 = abs(AveTotEnergy5_25t-AveTotEnergy5_0t)/abs(AveTotEnergy5_0t);
    AveErr10 = abs(AveTotEnergy10_25t-AveTotEnergy10_0t)/abs(AveTotEnergy10_0t);
    AveErr25 = abs(AveTotEnergy25_25t-AveTotEnergy25_0t)/abs(AveTotEnergy25_0t);
    AveErr100 = abs(AveTotEnergy100_25t-AveTotEnergy100_0t)/abs(AveTotEnergy100_0t);
    AveErr200 = abs(AveTotEnergy200_25t-AveTotEnergy200_0t)/abs(AveTotEnergy200_0t);

    figure; hold on; grid on;
    plot(1, AveErr1*100,'.','MarkerSize',25, 'DisplayName', '1 1/s'); 
    plot(5, AveErr5*100,'.','MarkerSize',25, 'DisplayName', '5 1/s'); 
    plot(10, AveErr10*100,'.','MarkerSize',25, 'DisplayName', '10 1/s'); 
    plot(25, AveErr25*100,'.','MarkerSize',25, 'DisplayName', '25 1/s'); 
    plot(100, AveErr100*100,'.','MarkerSize',25, 'DisplayName', '100 1/s'); 
    plot(200, AveErr200*100,'.','MarkerSize',25, 'DisplayName', '200 1/s');     
    xlabel("Strain Rate (1/s)");
    ylabel("Difference (Error) in Energy (%)");
    title("Difference in Internal Energy vs Applied Strain Rate");
    legend('show');

    figure; hold on; grid on;
    plot(1, AveTotEnergy_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(5, AveTotEnergy5_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(10, AveTotEnergy10_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(25, AveTotEnergy25_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(100, AveTotEnergy100_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    plot(200, AveTotEnergy200_0t,'b.','MarkerSize',25, 'DisplayName', '0 Fibers'); 
    
    plot(1, AveTotEnergy1_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(5, AveTotEnergy5_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(10, AveTotEnergy10_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(25, AveTotEnergy25_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers');
    plot(100, AveTotEnergy100_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    plot(200, AveTotEnergy200_25t,'r.','MarkerSize',25, 'DisplayName', '25 Fibers'); 
    xlabel("Strain Rate (1/s)");
    ylabel("Energy (J)");
    title("Internal Energy vs Applied Strain Rate");
    legend('show');
    
    
%% Vol Fraction
set(0,'defaultfigurecolor',[1 1 1]);

Abq1h_OG = ReadHost2("Ab_1h_Rate5");
[AbHost_1h_2t, AbTruss_1h_2t, AbqE_1h_2t] = ReadHostTruss2('Ab_1h_2t');
[AbHost_1h_25t, AbTruss_1h_25t, AbqE_1h_25t] = ReadHostTruss2('Ab_1h_25t');
[AbHost_1h_10t, AbTruss_1h_10t, AbqE_1h_10t] = ReadHostTruss2('Ab_1h_10t');

name1 = "1h 0t";
name2 = "1h 2t";
name3 = "1h 25t";
name4 = "1h 10t";

file1="explicit_3D";
name1 = "Flagshyp No Truss";

file2="embed_1h_2t";
file3="embed_1h_2t_correct";

FLAG_0 = ReadFlagshypOutputFile(file1, 83,1); 
FLAG_2 = ReadFlagshypOutputFile(file2, 82,1);
FLAG_2c = ReadFlagshypOutputFile(file3, 83,1);

file25="embed_1h_25t";
file25c="embed_1h_25t_correct";
 
FLAG_25 = ReadFlagshypOutputFile(file25, 82,1);
FLAG_25c = ReadFlagshypOutputFile(file25c, 82,1);

%% Kinetic Energy
    %Energy Difference (Error)    
    AveTotEnergy_0t = Abq1h_OG.KE(end);
    AveTotEnergy_2t = AbqE_1h_2t.KE(end);
    AveTotEnergy_25t = AbqE_1h_25t.KE(end);
    AveTotEnergy_10t = AbqE_1h_10t.KE(end);
    
    AveErr1 = abs(AveTotEnergy_0t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr2 = abs(AveTotEnergy_2t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr25 = abs(AveTotEnergy_25t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr10 = abs(AveTotEnergy_10t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);

    FlagEnergy_0t = FLAG_0.KE(end);
    FlagEnergy_2t = FLAG_2.KE(end);
    FlagEnergy_2tc = FLAG_2c.KE(end);
    FlagEnergy_25t = FLAG_25.KE(end);
    FlagEnergy_25tc = FLAG_25c.KE(end);
    
    FErr0 = abs(FlagEnergy_0t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr2 = abs(FlagEnergy_2t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr2c = abs(FlagEnergy_2tc-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr25 = abs(FlagEnergy_25t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr25c = abs(FlagEnergy_25tc-FlagEnergy_0t)/abs(FlagEnergy_0t);
    
    VolFrac0 = 0;
    VolFrac2 = 2*(0.02*1)/1;
    VolFrac25 = 25*(0.02*1)/1;
    VolFrac10 = 10*(0.02*1)/1;
    
Color0 = '#00008B';
Color1 = '#006400';
Color2 = '#B22222';
    
    figure; hold on; grid on;
    plot(VolFrac0, 100*AveErr1,'b.','MarkerSize',20, 'DisplayName', '0 Fibers-Abaqus'); 
        plot(VolFrac0, 100*FErr0,'b^','MarkerSize',10, 'DisplayName', '0 Fibers-Flagshyp'); 
    plot(VolFrac2, 100*AveErr2,'.','Color',Color2,'MarkerSize',20, 'DisplayName', '2 Fibers-Abaqus');
        plot(VolFrac2, 100*FErr2,'^','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp');
        plot(VolFrac2, 100*FErr2c,'x','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp Corrected');
%     plot(VolFrac10, 100*AveErr10,'^','MarkerSize',25, 'DisplayName', '10 Fibers-Abaqus');    
    plot(VolFrac25, 100*AveErr25,'.','Color',Color1,'MarkerSize',20, 'DisplayName', '25 Fibers-Abaqus'); 
        plot(VolFrac25, 100*FErr25,'^','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp');
        plot(VolFrac25, 100*FErr25c,'x','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp Corrected');
    xlabel("Volume Fraction");
    ylabel("Difference (Error) in Energy (%)");
    title("Difference in Kinetic Energy vs Fiber Volume Fraction");
    legend('show');

    figure; hold on; grid on;
    plot(VolFrac0, AveTotEnergy_0t,'b.','MarkerSize',20, 'DisplayName', '0 Fibers-Abaqus'); 
        plot(VolFrac0, FlagEnergy_0t,'b^','MarkerSize',10, 'DisplayName', '0 Fibers-Flagshyp'); 
    plot(VolFrac2, AveTotEnergy_2t,'.','Color',Color2,'MarkerSize',20, 'DisplayName', '2 Fibers-Abaqus');
        plot(VolFrac2, FlagEnergy_2t,'^','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp');
        plot(VolFrac2, FlagEnergy_2tc,'x','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp Corrected');
%     plot(VolFrac10, AveTotEnergy_10t,'^','MarkerSize',25, 'DisplayName', '10 Fibers-Abaqus');    
    plot(VolFrac25, AveTotEnergy_25t,'.','Color',Color1,'MarkerSize',20, 'DisplayName', '25 Fibers-Abaqus'); 
        plot(VolFrac25, FlagEnergy_25t,'^','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp');
        plot(VolFrac25, FlagEnergy_2tc,'x','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp Corrected');    
    xlabel("Volume Fraction");
    ylabel("Energy (J)");
    title("Kinetic Energy vs Fiber Volume Fraction");
    legend('show');

    
%% Internal Energy
    %Energy Difference (Error)    
    AveTotEnergy_0t = Abq1h_OG.IE(end);
    AveTotEnergy_2t = AbqE_1h_2t.IE(end);
    AveTotEnergy_25t = AbqE_1h_25t.IE(end);
    AveTotEnergy_10t = AbqE_1h_10t.IE(end);
    
    AveErr1 = abs(AveTotEnergy_0t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr2 = abs(AveTotEnergy_2t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr25 = abs(AveTotEnergy_25t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr10 = abs(AveTotEnergy_10t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    
    FlagEnergy_0t = FLAG_0.IE(end);
    FlagEnergy_2t = FLAG_2.IE(end);
    FlagEnergy_2tc = FLAG_2c.IE(end);
    FlagEnergy_25t = FLAG_25.IE(end);
    FlagEnergy_25tc = FLAG_25c.IE(end);
    
    FErr0 = abs(FlagEnergy_0t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr2 = abs(FlagEnergy_2t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr2c = abs(FlagEnergy_2tc-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr25 = abs(FlagEnergy_25t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr25c = abs(FlagEnergy_25tc-FlagEnergy_0t)/abs(FlagEnergy_0t);
    
    VolFrac0 = 0;
    VolFrac2 = 2*(0.02*1)/1;
    VolFrac25 = 25*(0.02*1)/1;
    VolFrac10 = 10*(0.02*1)/1;
    
    figure; hold on; grid on;
    plot(VolFrac0, 100*AveErr1,'b.','MarkerSize',20, 'DisplayName', '0 Fibers-Abaqus'); 
        plot(VolFrac0, 100*FErr0,'b^','MarkerSize',10, 'DisplayName', '0 Fibers-Flagshyp'); 
    plot(VolFrac2, 100*AveErr2,'.','Color',Color2,'MarkerSize',20, 'DisplayName', '2 Fibers-Abaqus');
        plot(VolFrac2, 100*FErr2,'^','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp');
        plot(VolFrac2, 100*FErr2c,'x','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp Corrected');
%     plot(VolFrac10, 100*AveErr10,'^','MarkerSize',25, 'DisplayName', '10 Fibers-Abaqus');    
    plot(VolFrac25, 100*AveErr25,'.','Color',Color1,'MarkerSize',20, 'DisplayName', '25 Fibers-Abaqus'); 
        plot(VolFrac25, 100*FErr25,'^','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp');
        plot(VolFrac25, 100*FErr25c,'x','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp Corrected');
    xlabel("Volume Fraction");
    ylabel("Difference (Error) in Energy (%)");
    title("Difference in Internal Energy vs Fiber Volume Fraction");
    legend('show');
    
    figure; hold on; grid on;
    plot(VolFrac0, AveTotEnergy_0t,'b.','MarkerSize',20, 'DisplayName', '0 Fibers-Abaqus'); 
        plot(VolFrac0, FlagEnergy_0t,'b^','MarkerSize',10, 'DisplayName', '0 Fibers-Flagshyp'); 
    plot(VolFrac2, AveTotEnergy_2t,'.','Color',Color2,'MarkerSize',20, 'DisplayName', '2 Fibers-Abaqus');
        plot(VolFrac2, FlagEnergy_2t,'^','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp');
        plot(VolFrac2, FlagEnergy_2tc,'x','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp Corrected');
%     plot(VolFrac10, AveTotEnergy_10t,'^','MarkerSize',25, 'DisplayName', '10 Fibers-Abaqus');    
    plot(VolFrac25, AveTotEnergy_25t,'.','Color',Color1,'MarkerSize',20, 'DisplayName', '25 Fibers-Abaqus'); 
        plot(VolFrac25, FlagEnergy_25t,'^','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp');
        plot(VolFrac25, FlagEnergy_2tc,'x','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp Corrected');   
    xlabel("Volume Fraction");
    ylabel("Energy (J)");
    title("Internal Energy vs Fiber Volume Fraction");
    legend('show');
    
%% Total Energy
    %Energy Difference (Error)    
    AveTotEnergy_0t = mean(Abq1h_OG.ETOTAL);
    AveTotEnergy_2t = mean(AbqE_1h_2t.ETOTAL);
    AveTotEnergy_25t = mean(AbqE_1h_25t.ETOTAL);
    AveTotEnergy_10t = mean(AbqE_1h_10t.ETOTAL);
    
    AveErr1 = abs(AveTotEnergy_0t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr2 = abs(AveTotEnergy_2t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr25 = abs(AveTotEnergy_25t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    AveErr10 = abs(AveTotEnergy_10t-AveTotEnergy_0t)/abs(AveTotEnergy_0t);
    
    FlagEnergy_0t = mean(FLAG_0.ET);
    FlagEnergy_2t = mean(FLAG_2.ET);
    FlagEnergy_2tc = mean(FLAG_2c.ET);
    FlagEnergy_25t = mean(FLAG_25.ET);
    FlagEnergy_25tc = mean(FLAG_25c.ET);
    
    FErr0 = abs(FlagEnergy_0t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr2 = abs(FlagEnergy_2t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr2c = abs(FlagEnergy_2tc-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr25 = abs(FlagEnergy_25t-FlagEnergy_0t)/abs(FlagEnergy_0t);
    FErr25c = abs(FlagEnergy_25tc-FlagEnergy_0t)/abs(FlagEnergy_0t);
    
    VolFrac0 = 0;
    VolFrac2 = 2*(0.02*1)/1;
    VolFrac25 = 25*(0.02*1)/1;
    VolFrac10 = 10*(0.02*1)/1;
    
    figure; hold on; grid on;
    plot(VolFrac0, 100*AveErr1,'b.','MarkerSize',20, 'DisplayName', '0 Fibers-Abaqus'); 
        plot(VolFrac0, 100*FErr0,'b^','MarkerSize',10, 'DisplayName', '0 Fibers-Flagshyp'); 
    plot(VolFrac2, 100*AveErr2,'.','Color',Color2,'MarkerSize',20, 'DisplayName', '2 Fibers-Abaqus');
        plot(VolFrac2, 100*FErr2,'^','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp');
        plot(VolFrac2, 100*FErr2c,'x','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp Corrected');
%     plot(VolFrac10, 100*AveErr10,'^','MarkerSize',25, 'DisplayName', '10 Fibers-Abaqus');    
    plot(VolFrac25, 100*AveErr25,'.','Color',Color1,'MarkerSize',20, 'DisplayName', '25 Fibers-Abaqus'); 
        plot(VolFrac25, 100*FErr25,'^','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp');
        plot(VolFrac25, 100*FErr25c,'x','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp Corrected');
    xlabel("Volume Fraction");
    ylabel("Difference (Error) in Energy (%)");
    title("Difference in Total Energy vs Fiber Volume Fraction");
    legend('show');
    
    figure; hold on; grid on;
    plot(VolFrac0, AveTotEnergy_0t,'b.','MarkerSize',20, 'DisplayName', '0 Fibers-Abaqus'); 
        plot(VolFrac0, FlagEnergy_0t,'b^','MarkerSize',10, 'DisplayName', '0 Fibers-Flagshyp'); 
    plot(VolFrac2, AveTotEnergy_2t,'.','Color',Color2,'MarkerSize',20, 'DisplayName', '2 Fibers-Abaqus');
        plot(VolFrac2, FlagEnergy_2t,'^','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp');
        plot(VolFrac2, FlagEnergy_2tc,'x','Color',Color2,'MarkerSize',10, 'DisplayName', '2 Fibers-Flagshyp Corrected');
%     plot(VolFrac10, AveTotEnergy_10t,'^','MarkerSize',25, 'DisplayName', '10 Fibers-Abaqus');    
    plot(VolFrac25, AveTotEnergy_25t,'.','Color',Color1,'MarkerSize',20, 'DisplayName', '25 Fibers-Abaqus'); 
        plot(VolFrac25, FlagEnergy_25t,'^','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp');
        plot(VolFrac25, FlagEnergy_2tc,'x','Color',Color1,'MarkerSize',10, 'DisplayName', '25 Fibers-Flagshyp Corrected');   
    xlabel("Volume Fraction");
    ylabel("Energy (J)");
    title("Total Energy vs Fiber Volume Fraction");
    legend('show');
%%
% Abq1h_OG = ReadHost2("Ab_1h_Force");
% [AbHost_1h_2t, AbTruss_1h_2t, AbqE_1h_2t] = ReadHostTruss2('Ab_1h_2t_Force');
% [AbHost_1h_25t, AbTruss_1h_25t, AbqE_1h_25t] = ReadHostTruss2('Ab_1h_25t_Force');
% 
% name1 = "1h 0t";
% name2 = "1h 2t";
% name3 = "1h 25t";
% 
% PlotEnergy3([Abq1h_OG.time, Abq1h_OG.KE], [AbqE_1h_2t.time, AbqE_1h_2t.KE], [AbqE_1h_25t.time, AbqE_1h_25t.KE], name1, name2,name3,'Kinetic Energy')
% PlotEnergy3([Abq1h_OG.time, Abq1h_OG.IE], [AbqE_1h_2t.time, AbqE_1h_2t.IE], [AbqE_1h_25t.time, AbqE_1h_25t.IE], name1, name2,name3,'Internal Energy')
% PlotEnergy3([Abq1h_OG.time, Abq1h_OG.WK], [AbqE_1h_2t.time, AbqE_1h_2t.WK], [AbqE_1h_25t.time, AbqE_1h_25t.WK], name1, name2,name3,'External Energy')
% PlotEnergy3([Abq1h_OG.time, Abq1h_OG.ETOTAL], [AbqE_1h_2t.time, AbqE_1h_2t.ETOTAL], [AbqE_1h_25t.time, AbqE_1h_25t.ETOTAL], name1, name2,name3,'Total Energy')
% 


%% Function Defs

function PlotEnergy(Data1, Data2, Name1, Name2,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy3(Data1, Data2, Data3, Name1, Name2,Name3,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    plot(Data3(:,1), Data3(:,2),'DisplayName',Name3,'LineWidth',1)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy5(Data1, Data2, Data3, Data4,Data5, Name1, Name2,Name3,Name4,Name5,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    plot(Data3(:,1), Data3(:,2),'DisplayName',Name3,'LineWidth',2)
    plot(Data4(:,1), Data4(:,2),'DisplayName',Name4,'LineWidth',1)
    plot(Data5(:,1), Data5(:,2),'DisplayName',Name5,'LineWidth',1)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy2(Data1, Data2, Name1, Name2,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end
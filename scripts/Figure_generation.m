% convergence_AP 
figure(1)
y_index = 0:1:30;
x_AP_3 = convergence_AP_list(1,:);
x_AP_4 = convergence_AP_list(2,:);
x_AP_5 = convergence_AP_list(3,:);
x_AP_6 = convergence_AP_list(4,:);

plot(y_index,x_AP_3,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_AP_4,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_AP_5,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_AP_6,"-x",'Markersize',10,'linewidth',3);

xlabel('Iteration index'),ylabel('Total energy consumption (J)');
legend('Proposed, M=3','Proposed, M=4', "Proposed, M=5", "Proposed, M=6");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

%convergence_q
figure(2)
y_index = 0:1:50;
x_q_10 = convergence_q_list(1,:);
x_q_50 = convergence_q_list(2,:);
x_q_100 = convergence_q_list(3,:);
x_q_200 = convergence_q_list(4,:);

plot(y_index,x_q_10,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_q_50,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_q_100,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_q_200,"-x",'Markersize',10,'linewidth',3);

xlabel('Iteration index'),ylabel('Total energy consumption (J)');
legend('Proposed, q=10','Proposed, q=50', "Proposed, q=100", "Proposed, q=200");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;


% Running time comparison

bar((run_time_list)')
x_value = [5,6,7,8,9,10];
set(gca,'XTickLabel',x_value);

xlabel('Number of users'),ylabel('Runing time (second)');
legend('Proposed, M=5','Proposed, M=10', "BCDO, AP=5", "BCDO, AP=10");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Monte Carlo offloading energy
figure(4)
y_index = 5:1:14;
x_OO = offloading_energy_list(1,:);
x_NO = offloading_energy_list(2,:);
x_RO = offloading_energy_list(3,:);
x_proposed = offloading_energy_list(4,:);
x_BCD = offloading_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);

xlabel('Number of users'),ylabel('Offloading energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Monte Carlo consensus energy
figure(4)
y_index = 5:1:14;
x_OO = consensus_energy_list(1,:);
x_NO = consensus_energy_list(2,:);
x_RO = consensus_energy_list(3,:);
x_proposed = consensus_energy_list(4,:);
x_BCD = consensus_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);

xlabel('Number of users'),ylabel('Consensus energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;


% Monte Carlo total energy
figure(5)
MC_energy_list = consensus_energy_list + offloading_energy_list;
y_index = 5:1:14;
x_OO = MC_energy_list(1,:);
x_NO =offloading_energy_list(2,:)  + consensus_energy_list(2,:) ;
x_RO = MC_energy_list(3,:);
x_proposed = MC_energy_list(4,:);
x_BCD = MC_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);

xlabel('Number of users'),ylabel('Total energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;




% Monte Carlo offloading delay
figure(6)
y_index = 5:1:14;
x_OO = offloading_delay_list(1,:);
x_NO = offloading_delay_list(2,:);
x_RO = offloading_delay_list(3,:);
x_proposed = offloading_delay_list(4,:);
x_BCD = offloading_delay_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);

xlabel('Number of users'),ylabel('Average offloading delay (ms)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Monte Carlo consensus delay 
y_index = 5:1:14;
x_OO = consensus_delay_list(1,:);
x_NO = consensus_delay_list(2,:);
x_RO = consensus_delay_list(3,:);
x_proposed = consensus_delay_list(4,:);
x_BCD = consensus_delay_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);

xlabel('Number of users'),ylabel('Consensus delay (ms)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Monte Carlo total delay 
MC_delay_list = consensus_delay_list + offloading_delay_list;
y_index = 5:1:14;
x_OO = MC_delay_list(1,:);
x_NO = MC_delay_list(2,:);
x_RO = MC_delay_list(3,:);
x_proposed = MC_delay_list(4,:);
x_BCD = MC_delay_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);

xlabel('Number of users'),ylabel('Total delay (ms)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Impact Paramaters of C 
bar((IP_C_energy_list)')
C_value = [10,20,30,40,50];
set(gca,'XTickLabel',C_value);
xlabel('Computing resources of APs (GHz)'),ylabel('Total energy consumption (J)');

legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;



% Impact Paramaters of B 
bar((IP_B_energy_list)')
B_value = [10,20,30,40,50];
set(gca,'XTickLabel',B_value);
xlabel('Bandwith of APs (MHz)'),ylabel('Total energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Impact Parameters of Sb
Sb_value = [100,200,300,400,500,600,700,800,900,1000];
y_index = Sb_value;
x_OO = IP_Sb_energy_list(1,:);
x_NO = IP_Sb_energy_list(2,:);
x_RO = IP_Sb_energy_list(3,:);
x_proposed = IP_Sb_energy_list(4,:);
x_BCD = IP_Sb_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);


xlabel('Size of block (kbits)'),ylabel('Total energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Impact Parameters of epsilon
epsilon_value = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5];
y_index = epsilon_value;
x_OO = IP_epsilon_energy_list(1,:);
x_NO = IP_epsilon_energy_list(2,:);
x_RO = IP_epsilon_energy_list(3,:);
x_proposed = IP_epsilon_energy_list(4,:);
x_BCD = IP_epsilon_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);


xlabel('Ratio of \epsilon_c to \epsilon_p'),ylabel('Total energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;


% Impact Parameters of AP
AP_value = [3,4,5,6,7,8,9,10,11,12];
y_index = AP_value;
x_OO = IP_AP_energy_list(1,:);
x_NO = IP_AP_energy_list(2,:);
x_RO = IP_AP_energy_list(3,:);
x_proposed = IP_AP_energy_list(4,:);
x_BCD = IP_AP_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);


xlabel('Number of APs'),ylabel('Total energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;

% Impact Parameters of Computing density
Computing_value = [500,1000,1500,2000,2500,3000,3500,4000,4500,5000];
y_index = Computing_value;
x_OO = IP_computing_energy_list(1,:);
x_NO = IP_computing_energy_list(2,:);
x_RO = IP_computing_energy_list(3,:);
x_proposed = IP_computing_energy_list(4,:);
x_BCD = IP_computing_energy_list(5,:);

plot(y_index,x_OO,"-o",'Markersize',10,'linewidth',3); hold on;
plot(y_index,x_NO,"-d",'Markersize',10,'linewidth',3);
plot(y_index,x_RO,"-*",'Markersize',10,'linewidth',3);
plot(y_index,x_proposed,"-x",'Markersize',10,'linewidth',3);
plot(y_index,x_BCD,"-^",'Markersize',10,'linewidth',3);


xlabel('Computing density (CPU cycles/bit)'),ylabel('Total energy consumption (J)');
legend('OO','SO', "RO", "Proposed", "BCDO");
set(gca,'FontName','Times New Roman','FontSize',12);
grid on;



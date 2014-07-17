data = [0	0.6	1	1.4	1.6	1.8	2	2.2	2.4	2.6	3	3.4	4	0.84	1.12	1.52	1.66	2.12	2.68	2.96;...
14.9	16.8	34.8	141	152	287	904	990	531	118	162	5.41	6.09	71.2	19.4	170	151	1000	92	168];

dataSorted = sortrows(data',1)';
figure()
plot(dataSorted(1,:)*1.023,dataSorted(2,:)/1180, 'linewidth', 4);

xlabel('incidence angle / degrees')
ylabel('efficiency')
grid minor
set(gcf,'color', 'w')  
set(gca,'fontsize', 24) 

xlbl = get(gca,'xlabel')
set(xlbl,'fontsize', 30) 

ylbl = get(gca,'ylabel')
set(ylbl,'fontsize', 30) 

xlim([0 4])
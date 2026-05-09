%clear all
close all
clc

% eseguire prima Biellittica_bitangente per i dati
figure
plot(delta_v(end,:,1),delta_t(end,:,1),"k--",delta_v(end,:,2),delta_t(end,:,2),"k")
grid on
xlabel("Δv (km/s)",'fontsize',12)
ylabel("Δt (s)",'fontsize',12)
set(gcf,'color','w');
hold on

% manovra standard
STDx = [8.7215];
STDy = [31661.1];
STD = scatter(STDx,STDy, 30,'d', "red", "filled");

% variazione standard 1

STD1x = [7.7063];
STD1y = [41259.4];

STD1 = scatter(STD1x,STD1y, 30,'d', "blue", "filled");

% variazione standard 2
STD2x = [6.7227];
STD2y = [23747.36];

STD2 = scatter(STD2x,STD2y, 30,'d', "green", "filled");

% diretta
DRx = [18.9954];
DRy = [1980.6];

DIRECT = scatter(DRx,DRy, 30,'d', "black", "filled");

% 2-impulsi best trade off
DR1x = [5.8704];
DR1y = [15722.15];

DIRECT1 = scatter(DR1x,DR1y, 30,'d', "cyan", "filled");

% biellittica bitangente
BBx = [4.8466];
BBy = [153342.67];

BIELLIPTIC = scatter(BBx,BBy, 30,'d', "magenta", "filled");

grid on

legend('Bielliptic A','Bielliptic B','Std','Std1','Std2','Direct1','Direct2','Chosen bielliptic','fontsize',9,'Location','best')
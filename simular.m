% Corre las simulaciones con el modelo SIR
% para los diferentes escenarios de R0
% @ktakahashig
% 29 de abril de 2020
%

clear

dt=1; % Paso de tiempo (dias)
tzero=datenum(2020,1,1); % t = 0
tini=datenum(2020,3,4); % Inicio de brote
t0sim=datenum(2020,3,15); % Inicio de simulacion
I0=1046.6; % Numero de infectados a inicio de simulacion
t=[tini:tini+365]'; % Periodo total

i0=find(t==t0sim); % Auxiliares
Nt=length(t);

% Parametros de modelo
N=32.6e6; % Poblacion nacional
gamma = 0.111; % Estimado en la calibracion

% Cargar escenarios de R0 base y la estacionalidad
load R0base.txt; load R0est.txt; 
R0est=[zeros(Nt,1) R0est]; % Incluir estacionalidad cero
Nbase=size(R0base,2); Nest=size(R0est,2);

% Correr el modelo para combinaciones de escenario base y estacionalidad
SS=[]; II=[]; RR=[];
for m=1:Nbase
   for n=1:Nest
      I=I0; S=N-I; R=0; X=[S I R]; 
      for i=i0+1:Nt
             beta=(R0base(i,m)+R0est(i,n))*gamma; % beta para este tiempo
             X=[X;rk4(X(end,:),dt,[beta gamma N])]; % Paso de modelo
      end
      SS=[SS X(:,1)]; II=[II X(:,2)]; RR=[RR X(:,3)]; % Juntar salidas
   end
end

% Guardar salidas
save -ascii simS.txt SS
save -ascii simI.txt II
save -ascii simR.txt RR



% Figuras  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsim=[t0sim:t(end)]';
labs=strvcat('a','b','c','d');
scen=strvcat('A','B','C','D');
option=strvcat('k-','r-.','g:','c--','b--','m-');


figure(1,'papersize',[9 4.5]),clf
for m=1:Nbase
   subplot(1,Nbase,m)
   h=[];
   for i=1:Nest
      ii=Nest*(m-1)+i;
      h=[h;semilogy(tsim,II(:,ii),option(i,:))];
      hold on
   end
   hold off
   xlim([t(1) t(end)]);
   ylim([1e1 1e7])
   set(h(1),'linewidth',2.8)
   set(h(2:end),'linewidth',2.5);
   set(h(end),'linewidth',2);
   %set(gca,'XTick',datenum(2020,[3:12]',1))
   set(gca,'XTick',datenum([repmat(2020,10,1); repmat(2021,4,1)],[[3:12] [1:4]]',1))
   datetick('mm','keeplimits','keepticks')
   set(gca,'fontsize',10)
   xlabel('Mes','fontsize',11)
   ylabel('Numero de personas','fontsize',11)
   title([labs(m) '\) Infectados: Escenarios ' scen(m)],'fontsize',13)
end
legend('Base','Kissler','Neher','Wang 1','Wang 2','Merow&Urban','location','southwest')
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','simI.pdf')


figure(2,'papersize',[9 4.5]),clf
for m=1:Nbase
   subplot(1,Nbase,m)
   h=[];
   for i=1:Nest
      ii=Nest*(m-1)+i;
      h=[h;semilogy(tsim,II(:,ii)+RR(:,ii),option(i,:))];
      hold on
   end
   hold off
   xlim([t(1) t(end)]);
   ylim([1e3 1e8]);
   set(h(1),'linewidth',2.8)
   set(h(2:end),'linewidth',2.5);
   set(h(end),'linewidth',2);
   %set(gca,'XTick',datenum(2020,[3:12]',1))
   set(gca,'XTick',datenum([repmat(2020,10,1); repmat(2021,4,1)],[[3:12] [1:4]]',1))
   datetick('mm','keeplimits','keepticks')
   set(gca,'fontsize',9)
   xlabel('Mes','fontsize',10)
   ylabel('Numero de personas','fontsize',10)
   title([labs(m) '\) Casos totales (I+R): Escenarios ' scen(m)],'fontsize',13)
end
legend('Base','Kissler','Neher','Wang 1','Wang 2','Merow&Urban','location','southeast')
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','simIR.pdf')


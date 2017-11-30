%
close all;
clear;
clc;

%====================================%
% simulation data
%====================================%

skip=0;
nrSamples = 0;
Ainlet=1.000000e-02;    % cross sectional area in m2
Utarget=10;            % target superficial velocity in m/s
Vtarget=Utarget*Ainlet; % target flowrate in m3/s

filename = 'surfaceFieldValue.dat';
timedirectory = '0';

    %Daten einlesen
    path = '../../CFD/postProcessing/volFlow_inlet';
    data = transpose(load(strcat(path,'/',timedirectory,'/',filename)));          %data
    [x,y]=size(data);
    time=data(1,1+skip:y)
    volFlow_inlet=data(2,1+skip:y)

    %Daten einlesen
    path = '../../CFD/postProcessing/volFlow_outlet';
    data = transpose(load(strcat(path,'/',timedirectory,'/',filename)));          %data
    [x,y]=size(data)
    time=data(1,1+skip:y);
    volFlow_outlet=data(2,1+skip:y);

    %Daten einlesen
    path = '../../CFD/postProcessing/volFlow_wall';
    data = transpose(load(strcat(path,'/',timedirectory,'/',filename)));          %data
    [x,y]=size(data)
    time=data(1,1+skip:y);
    volFlow_wall=data(2,1+skip:y);

%- time integrate data to vol entered per TS
deltatT=time(2)-time(1)
volPerTs_inlet=volFlow_inlet(1,:).*deltatT;
volPerTs_outlet=volFlow_outlet(1,:).*deltatT;
volPerTs_wall=volFlow_wall(1,:).*deltatT;

%- accumulated vol entered
vol_inlet(1)=volPerTs_inlet(1);
vol_outlet(1)=volPerTs_outlet(1);
vol_wall(1)=volPerTs_wall(1);
for i=2:y-skip
    vol_inlet(i) = vol_inlet(i-1) + volPerTs_inlet(i);
    vol_outlet(i) = vol_outlet(i-1) + volPerTs_outlet(i);
    vol_wall(i) = vol_wall(i-1) + volPerTs_wall(i);
end

%===================================
% plot 1
xAxisLabel = 'time [s]';
yAxisLabel = 'time integrated flux in [m3]';

% Create figure
figure1 = figure('PaperPositionMode','manual','PaperUnits','centimeters',...
%    'PaperPosition',[0 0 15.5 10],'PaperSize',[15.5 10],...
    'Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','LineWidth',1,...
    'FontWeight','normal','FontSize',11,'FontName','Helvetica-Narrow');
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(time,vol_inlet,'r','Parent',axes1,'Marker','none','LineWidth',1,...
     time,-vol_outlet,'g','Parent',axes1,'Marker','none','LineWidth',1,...
     time,-vol_wall,'b','Parent',axes1,'Marker','none','LineWidth',1);

% Create xlabel
xlabel(xAxisLabel,'FontWeight','bold','FontSize',11,'FontName','Helvetica-Narrow');
% Create ylabel
ylabel(yAxisLabel,'FontWeight','bold','FontSize',11,'FontName','Helvetica-Narrow');
% Define axis
%axis([dMin,dMax]);

title('time integrated flux acoss patches');
legend('inlet','-1*outlet','wall');
print -color "volEntered.png"

%===================================
% plot 2
xAxisLabel = 'time [s]';
yAxisLabel = 'vol entered in TS [m3]';

% Create figure
figure2 = figure('PaperPositionMode','manual','PaperUnits','centimeters',...
%    'PaperPosition',[0 0 15.5 10],'PaperSize',[15.5 10],...
    'Color',[1 1 1]);

% Create axes
axes2 = axes('Parent',figure2,'YGrid','on','XGrid','on','LineWidth',1,...
    'FontWeight','normal','FontSize',11,'FontName','Helvetica-Narrow');
box(axes2,'on');
hold(axes2,'all');

% Create plot
plot(time,volPerTs_inlet,'r','Parent',axes2,'Marker','none','LineWidth',1,...
     time,-volPerTs_outlet,'g','Parent',axes2,'Marker','none','LineWidth',1,...
     time,-volPerTs_wall,'b','Parent',axes2,'Marker','none','LineWidth',1);

% Create xlabel
xlabel(xAxisLabel,'FontWeight','bold','FontSize',11,'FontName','Helvetica-Narrow');
% Create ylabel
ylabel(yAxisLabel,'FontWeight','bold','FontSize',11,'FontName','Helvetica-Narrow');
% Define axis
%axis([dMin,dMax]);

title('vol flux * time = vol entered per TS');
legend('inlet','outlet','wall');
print -color "volEnteredPerTS.png"

%===================================
% plot 3
xAxisLabel = 'time [s]';
yAxisLabel = 'vol flux in [m3/s]';

% Create figure
figure3 = figure('PaperPositionMode','manual','PaperUnits','centimeters',...
%    'PaperPosition',[0 0 15.5 10],'PaperSize',[15.5 10],...
    'Color',[1 1 1]);

% Create axes
axes3 = axes('Parent',figure3,'YGrid','on','XGrid','on','LineWidth',1,...
    'FontWeight','normal','FontSize',11,'FontName','Helvetica-Narrow');
box(axes3,'on');
hold(axes3,'all');

% Create plot
plot(time,volFlow_inlet(1,:),'r','Parent',axes3,'Marker','none','LineWidth',1,...
     time,-volFlow_outlet(1,:),'g','Parent',axes3,'Marker','none','LineWidth',1,...
     time,-volFlow_wall(1,:),'b','Parent',axes3,'Marker','none','LineWidth',1,...
     [min(time),max(time)],[-Vtarget,-Vtarget],'c','Parent',axes3,'Marker','none','LineWidth',1);

% Create xlabel
xlabel(xAxisLabel,'FontWeight','bold','FontSize',11,'FontName','Helvetica-Narrow');
% Create ylabel
ylabel(yAxisLabel,'FontWeight','bold','FontSize',11,'FontName','Helvetica-Narrow');
% Define axis
%axis([dMin,dMax]);

title('vol flux');
legend('inlet','outlet','wall','Vtarget');
print -color "volflow.png"


% Script to plot distance vs time seismograms

stations_lon = {'LON2', 'LON4', 'LON6', 'LON8', 'LON10'};
stations_lat = {'LAT0', 'LAT1', 'LAT2', 'LAT3', 'LAT4'}';

%% Plot stacked
for i=(1:length(stations_lat))
    afile = ['raw_data/sim1/II.',char(stations_lat(i)),'.RTZ.ascii'];
%     afile2 = ['raw_data/sim3/II.',char(station2(i)),'.RTZ.ascii'];
    
    adata = load(afile);
%     adata2 = load(afile2);
    
    ta = adata(:,1);
    
    da = adata(:,4);
%     da2 = adata2(:,4);
    
    if i==1
        figure('Position',[450 360 1100 470])
        plot(ta,da,'LineWidth',1)
        dx = 1.0*max(da);
        
        xlabel('time [s]')
        ylabel('distance from source in x-direction [km]')
        ylim([-1.2*max(da) (length(stations_lat)-0.3)*dx])
        
        hold on
        %plot(ta+0.25,da,'color',[0.301 0.745 0.933])
        %plot(ts,ds2,'r--')
        
    else
        plot(ta,da + (i-1)*dx,'LineWidth',1)
        %plot(ta+0.25,da + (i-1)*dx,'color',[0.301 0.745 0.933])
        %plot(ta+0.25,da2 + (i-1)*dx,'r--')
        %plot(ts,ds2 + (i-1)*dx,'r--')
    end

end

set(gca,'YTick',(0:dx:(length(stations_lat)-1)*dx),'YTickLabel',{'0','22','44','66','88'}')
xlim([0.25 4])
set(gca, 'YDir')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Plot just one 
figure(10)
plot(ta,da,'k','LineWidth',1)
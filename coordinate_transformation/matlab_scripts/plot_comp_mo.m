% Script to plot distance vs time seismograms

station = {'00','10','20','30','40'};
station2 = {'00','01','02','03','04'};

for i=(1:length(station))
    afile = ['II.station',char(station(i)),'.RTZ.ascii'];
    afile2 = ['II.station',char(station2(i)),'.RTZ.ascii'];
    
    sdata = dlmread(sfile,',',1,0);
    sdata2 = dlmread(sfile2,',',1,0);
    adata = load(afile);
    adata2 = load(afile2);
    
    ts = [0; sdata(:,1)];
    ta = adata(:,1);
    
    ds = [0; sdata(:,2)];
    ds2 = [0; sdata2(:,2)];
    da = adata(:,4);
    da2 = adata2(:,4);
    
    if i==1
        figure('Position',[450 360 1100 470])
        plot(ts,ds,'k','LineWidth',1)
        dx = 0.8*max(ds);
        
        xlabel('time [s]')
        ylabel('distance from source in x-direction [km]')
        ylim([-1.2*max(ds) (length(station)-0.3)*dx])
        
        hold on
        plot(ta+0.25,da,'r-.')
        %plot(ta+0.25,da,'color',[0.301 0.745 0.933])
        %plot(ts,ds2,'r--')
        
    else
        plot(ts,ds + (i-1)*dx,'k','LineWidth',1)
        plot(ta+0.25,da + (i-1)*dx,'r-.')
        %plot(ta+0.25,da + (i-1)*dx,'color',[0.301 0.745 0.933])
        %plot(ta+0.25,da2 + (i-1)*dx,'r--')
        %plot(ts,ds2 + (i-1)*dx,'r--')
    end

end

set(gca,'YTick',(0:dx:(length(station)-1)*dx),'YTickLabel',{'0.0','0.5','1.0','1.5','2.0'})
xlim([0.25 4])
set(gca, 'YDir', 'reverse')
legend('SPECFEM3D', 'AxiSEM3D')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
figure
times=load('raw_data/geographic_bathymetry/wavefields/times.txt');
times=round(unique(times),1);

for f_ = (0:length(times)-1)
    f=f_;%*floor(length(times)/10);
    data = load(['raw_data/geographic_bathymetry/wavefields/Phi_0.000000_Frame',int2str(f),'.txt']); % x z u_s u_p u_z
    data = sortrows(data);
    
    x = data(:,1);
    z = data(:,2) + data(:,3);
    u = data(:,6);
    
    subplot(2,5,f_+1)
    scatter(x,z,200,u,'.')
    
    title(['t = ', num2str(times(f+1)), 's'])
    caxis([-5 5] * 10^(-4))
end
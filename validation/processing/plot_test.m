%% Test the outputs of simulation 1 
% Cartesian equivalent of domain as calculated in 
% get_dom_sizes.py

data1_lat2 = load('./raw_data/sim1/II.LAT2.RTZ.ascii');
data1_lat3 = load('./raw_data/sim1/II.LAT3.RTZ.ascii');
data1_lat4 = load('./raw_data/sim1/II.LAT4.RTZ.ascii');

figure(1)
subplot(3,1,1)
plot(data1_lat2(1:100,1), data1_lat2(1:100,4), 'b-')
xlabel('Time (sec)')
ylabel('Radial displacement (mm)')
title('At source')

subplot(3,1,2)
plot(data1_lat3(1:500,1), data1_lat3(1:500,4), 'b-')
xlabel('Time (sec)')
ylabel('Radial displacement (mm)')
title('20 km away')

subplot(3,1,3)
plot(data1_lat4(1:1000,1), data1_lat4(1:1000,4), 'b-')
xlabel('Time (sec)')
ylabel('Radial displacement (mm)')
title('30 km away')
%% Test the outputs of simulations against one another
% 1. Cartesian equivalent (calculated assuming spherical Earth and flat
% domain)
% 2. Spherical equivalent (calculated assuming spherical Earth defined by
% the radius in the centre of the domain)
% 3. Geographic domain (obtained by transforming into geocentric
% coordinates, and then calculated by calculating the length of lat and lon
% degree at each lat)

data1_lat9 = load('./raw_data/sim1/II.LAT4.RTZ.ascii');
data2_lat9 = load('./raw_data/sim2/II.LAT4.RTZ.ascii');
data3_lat9 = load('./raw_data/sim3/II.LAT4.RTZ.ascii');

figure(2)
subplot(3,1,1)
plot(data1_lat9(300:500,1), data1_lat9(300:500,4), 'b-')
xlabel('Time (sec)')
title('Cartesian')

subplot(3,1,2)
plot(data2_lat9(300:500,1), data2_lat9(300:500,4), 'b-')
xlabel('Time (sec)')
title('Spherical')

subplot(3,1,3)
plot(data3_lat9(300:500,1), data3_lat9(300:500,4), 'b-')
xlabel('Time (sec)')
title('Geographic')
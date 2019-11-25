% queryCRUST1.m   This script should be run to query the model CRUST1.0
% at specific locations (one point at a time) by prompting for latitude
% (lat) and longitude (lon). This script needs access to the functions and 
% the datastructure found in the CRUST1 toolbox. Obviously that tolbox 
% needs to be on the user's path, which can arranged using matlab's
% native functions pathtool.m or addpath.m

% Version 1.1            Michael Bevis          21 June 2017

clear all
fprintf(1,'loading CRUST1.0 structure C1\n')
load('CRUST1.mat')
if ~exist('C1')
    error('Could not load CRUST1.mat')
end

% query user for lat,lon, identify the cell(s) and report the profile(s)
fin=0;
while fin==0
  accept1=0;
  while accept1==0
        lat=input('Enter lat  (or just return to quit) : ');
        if isempty(lat)
            return
        elseif lat>= -90 & lat<=90       
            accept1=1;
        end
  end
  accept2=0;
  while accept2==0 & ~isempty(lat)
        lon=input('Enter lon  (or just return to quit) : ');
        if isempty(lon)
            return
        elseif lon>=-180 & lon<=180
            accept2=1;
        else
            accept2=0;
        end
  end
  if ~isempty(lat) && ~isempty(lon) 
      fprintf(1,'\nQuery CRUST1.0 at lat=%f   lon=%f  \n',[lat lon]);
      [i,j]= findCrust1cell(lat,lon);
      ncell=length(i);
      if ncell>1
            fprintf(1,'This pt is associated with %i CRUST1.0 cells\n',...
                      ncell);
      end
      for k=1:ncell
        clat=C1.CLAT(i(k),j(k));
        clon=C1.CLON(i(k),j(k));
        fprintf(1,'\nCRUST1.0 CELL CENTERED AT CLAT = %5.1f ',clat)
        fprintf(1,' CLON = %6.1f\n',clon)
        [b,p,s,r] = getCrust1(i(k),j(k),C1);
        %block=[b p s r]
        t=b(1:end-1)-b(2:end); % layer thickness (layers 1-8
        topo=b(1);
        fprintf(1,'mean surface height: %9.3f km\n',topo)
        db=diff(b);
        kk=find(db~=0);
        nlayer=length(kk);
        fprintf(1,['-- LAYER TYPE ---       Vp    Vs    rho    bottom',...
                  '    thick     mu     lambda     E       nu \n'])
        fprintf(1,['                       km/s  km/s  g/cm3     km',...
                   '       km      GPa      GPa      Gpa \n']);
        for ii=1:nlayer 
            n=kk(ii);   % layer number
            fprintf(1,'%17s  ',layer{n})
            fprintf(1,'    %4.2f ',p(n))
            fprintf(1,' %4.2f ',s(n))
            fprintf(1,' %5.2f  ',r(n))
            fprintf(1,'%7.2f ',b(n+1))
            fprintf(1,'  %6.2f ',t(n))
            rho=r(n)*1e3; % SI units kg/m3
            Vp=p(n)*1e3;  % SI units m/s
            Vs=s(n)*1e3;  % SI units m/s
            mu=rho*Vs^2;  % SI units Pa
            lambda=rho*(Vp^2-2*Vs^2);  % SI units Pa
            E = (mu.*(3*lambda+2*mu))./(lambda+mu);
            nu= lambda./(2*(lambda+mu));
            fprintf(1,' %7.3f ',mu*1e-9)
            fprintf(1,' %7.3f ',lambda*1e-9)
            fprintf(1,' %7.3f ',E*1e-9)
            fprintf(1,' %5.3f ',nu) 
            fprintf(1,'\n')
        end
        fprintf(1,'%17s  ','uppermost mantle')
        fprintf(1,'    %4.2f ',p(end))
        fprintf(1,' %4.2f ',s(end))
        fprintf(1,' %5.2f  ',r(end))
        fprintf(1,'%s',blanks(17))
        rho=r(end)*1e3; % SI units kg/m3
        Vp=p(end)*1e3;  % SI units m/s
        Vs=s(end)*1e3;  % SI units m/s
        mu=rho*Vs^2;  % SI units Pa
        lambda=rho*(Vp^2-2*Vs^2);  % SI units Pa
        E = (mu.*(3*lambda+2*mu))./(lambda+mu);
        nu= lambda./(2*(lambda+mu));
        fprintf(1,' %7.3f ',mu*1e-9)
        fprintf(1,' %7.3f ',lambda*1e-9)
        fprintf(1,' %7.3f ',E*1e-9)
        fprintf(1,' %5.3f ',nu) 
        fprintf(1,'\n')
      end
      fprintf(1,'\n\n')
  else
      fin=1;
  end
end
    
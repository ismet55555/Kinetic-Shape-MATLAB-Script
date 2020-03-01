clc         %Clear the command window
close all   %Close all open figures/matlab windows
clear all   %Clear all data in matlab workspace

%==============================================================================
%                                    INPUT 
%==============================================================================
n          = 5;  %Numerical Resolution
Resolution = 100;

%Definition Range (Radians)
theta1 = linspace(0, 2*pi ,Resolution*n);%Elevation Angle Range (theta)
phi1 = linspace(0, 2*pi ,Resolution*n);  %Azimuth Angle Range (phi)
[theta, phi] = meshgrid(theta1,phi1);

%Applied vertical force to axle of 3D KS  (Unit Force)
%  Note:  This is a constant funciton. While this code does not support a 
%         variable vertical force funciton, it can be added.
Fv = 1000;

%Generating Desired Reaction Force Function in Theta (raidal) Direction 
%  Note: This code only supports horizontal ground reaction forces in 
%        exponential and sinosodial force functions, specifiying constants.
c1R = 100;      odrR = 0;
c2R = 5;        odr2R = 1; 
c3R = 10;       freq1R = 1;
c4R = 0;        freq2R = 1;

%Fr = zeros(length(theta)); %Blank Vector
Fr =   c1R*theta.^odrR       ...
     + c2R*phi.^odr2R        ...
     + c3R*sin(freq1R*theta) ...
     + c4R*sin(freq2R*phi);

%Generating Desired Reaction Force Function in Phi (tangential) Direction 
c1T = 10;      odrT = 0;
c2T = 10;      odr2T = 0;
c3T = 0;       freq1T = 1;
c4T = 10;      freq2T = 1;

%FT = zeros(length(phi)); %Blank Vector
FT =  c1T*theta.^odrT        ...
    + c2T*phi.^odr2T        ...
    + c3T*sin(freq1T*theta) ...
    + c4T*sin(freq2T*phi);

%Specify coordinates where within Input force the shape will start
%     Coordinates = [Theta, Phi]
pThidx = 1*n;   %Starting Point Index in theta direction
pPhidx = 1*n;   %Starting Point Index in phi direction

%Initial Radius of 3D Kinetic Shape  (Unit Length)
Ri = 1;

%Integration direction values
%     A value of 1 is 2pi, 0.5 is pi, and so on...
ThetaPath = 1; 
PhiPath   = 1;
%==============================================================================


FTotal = Fr + FT;       %Total Force - Radial plus Tangential Force

%Direction of Integration (Rolling Movement)
%     In this code the integration direciton is linear and constantly
%     increasing, however one can reprogram this code such that the
%     integration path is over a non-linear theta and phi range.
position = [[pThidx, pPhidx], [theta1(pThidx), phi1(pPhidx)], FTotal(pPhidx, pThidx)];
for r = 8:(Resolution*n)   
    %Direction travel 
    %For each angle step, take one index step in theta direction
    pThidx = pThidx + ceil(ThetaPath*4);   %Index amount moved in theta direction
    pPhidx = pPhidx + ceil(PhiPath*4);   %Index amount moved in phi direction
    
    pThidx = ceil(r*ThetaPath);
    pPhidx = ceil(r*PhiPath);

    position = [position; pThidx pPhidx theta1(pThidx), phi1(pPhidx), FTotal(pPhidx, pThidx)];
end



%==============================================================================
%               Applying Three-Dimensional Kinetic Shape Equation 
%==============================================================================
%Derivation of Shape (Curve)
[a, b] = size(theta);
r1(a,b) = 0;
firstRun = 0;
for j = 2:length(position)
    t      = position(j,   1);
    t_min1 = position(j-1, 1);
    
    p      = position(j,   2);
    p_min1 = position(j-1, 2);

    th  = phi(t,      p);
    th1 = phi(t_min1, p);
    
    ph  = theta(t, p);
    ph1 = theta(t, p_min1);
    
    
    if(firstRun==0)
        r1(t_min1,p_min1) =  Ri;
        firstRun = 1;
    end
    
    %==============  Theta Direction Integrals =================
    Int_Fr_dth = c1R*(th^(odrR+1))/(odrR+1) ... 
       +      th*c2R*(ph^(odr2R))...
               - c3R*(1/freq1R)*cos(freq1R*th) ...
       +      th*c4R*sin(freq1R*ph);
    
    Int_Fr_dth1 = c1R*(th1^(odrR+1))/(odrR+1)  ...
           +  th1*c2R*(ph^(odr2R))...
                - c3R*(1/freq1R)*cos(freq1R*th1) ...
           +  th1*c4R*sin(freq1R*ph);
    C = log( r1(t_min1,p_min1) ); %Constant

    r1(t_min1,p) = exp(  (Int_Fr_dth / Fv) - (Int_Fr_dth1 / Fv)   +  C  );
    
    
    %==================  Phi Direction Integrals ===================
    Int_FT_dph =         ph*c1T*th^odrT + c2T*(ph^(odr2T+1))/(odr2T+1)...
               +  ph*c3T*sin(freq1T*th) - c4T*(1/freq2T)*cos(freq2T*ph);
    Int_FT_dph1 =       ph1*c1T*th^odrT + c2T*(ph1^(odr2T+1))/(odr2T+1)...
              +  ph1*c3T*cos(freq1T*th) - c4T*(1/freq2T)*cos(freq2T*ph1);
    
    Int_Fr_dph =         ph*c1R*th^odrR + c2R*(ph^(odr2R+1))/(odr2R+1)...
                + ph*c3R*sin(freq1R*th) - c4R*(1/freq2R)*cos(freq2R*ph);
    Int_Fr_dph1 =       ph*c1R*th1^odrR + c2R*(ph1^(odr2R+1))/(odr2R+1)...
               + ph1*c3R*sin(freq1R*th) - c4R*(1/freq2R)*cos(freq2R*ph1);
       
    r1(t,p) = r1(t_min1,p)* exp(  (Int_FT_dph/Int_Fr_dph) ...
              - (Int_FT_dph1/Int_Fr_dph1)    );
end

%==============================================================================
%                           Plotting Kinetic Shape 
%==============================================================================
%Converting Shape from spherical to cartesian coordiantes for plotting
[x1,y1,z1] = sph2cart(theta,phi,r1);
value = find(x1 ~= 0);
noZigZag = [x1(value(2:2:end)),y1(value(2:2:end)),z1(value(2:2:end))];

%Plotting Derived 3D KS
figure(1);
plot3(0,0,0,'-mo','MarkerSize',10,'MarkerFaceColor','k')   %Origin point
hold on
plot3(x1(value(2:2:end)),y1(value(2:2:end)),z1(value(2:2:end)),...
      '.k','MarkerSize',10);
title('3D Kinetic Shape', 'FontSize',14, 'FontName','Times New Roman')
xlabel('X') ; ylabel('Y') ; zlabel('Z')
axis equal %tight
view(60, 30);
grid on

%Plotting Radial Force
figure(2);  %Input
contour(theta1, phi1, Fr, 18);
hold on
plot(position(:,3), position(:,4),'--k','Linewidth',4)
title('Radial Reaction Force', 'FontSize',14, 'FontName','Times New Roman')
xlabel('Theta (Radians)','FontSize',12, 'FontName','Times New Roman'); 
ylabel('Phi (Radians)','FontSize',12, 'FontName','Times New Roman'), 
zlabel('Force','FontSize',12, 'FontName','Times New Roman');
set(gca,'xLim',[min(theta1),max(theta1)]);
set(gca,'yLim',[min(phi1),max(phi1)]);
colorbar
grid on

%Plotting Tangential Force - Input
figure(3);  
contour(theta1, phi1, FT, 18);
hold on
plot(position(:,3), position(:,4),'--k','Linewidth',4)
title('Tangential Reaction Force', 'FontSize',14, 'FontName','Times New Roman')
xlabel('Theta (Radians)','FontSize',12, 'FontName','Times New Roman'); 
ylabel('Phi (Radians)','FontSize',12, 'FontName','Times New Roman'), 
zlabel('Force','FontSize',12, 'FontName','Times New Roman');
set(gca,'xLim',[min(theta1),max(theta1)]);
set(gca,'yLim',[min(phi1),max(phi1)]);
colorbar
grid on

%Ploting Total Forces onto a contour plot
figure(4);
C = contour(theta1, phi1, FTotal, 18);
hold on
plot(position(:,3), position(:,4),'--k','Linewidth',4)
title('TOTAL Reaction Force Input', 'FontSize',14, 'FontName','Times New Roman')
xlabel('Theta (Radians)','FontSize',12, 'FontName','Times New Roman');  
ylabel('Phi (Radians)','FontSize',12, 'FontName','Times New Roman')
set(gca,'xLim',[min(theta1),max(theta1)]);
set(gca,'yLim',[min(phi1),max(phi1)]);
colorbar
grid on

%==============================================================================
%            Export 3D data points to text file and .STL File
%==============================================================================
%XYZ point coordinates to text file:
Data = noZigZag;
dlmwrite('3D_KS_XYZ.txt', Data, 'precision', '%g', ...
          'newline', 'pc');
      
%3D object to .stl file:
surf2stl('3D_KS_STL.stl',x1,y1,z1)
      
%Checking what was written to the .stl file file
check = stlread('3D_KS_STL.stl')
figure
patch(check,'FaceColor',       [0.5 0.8 1.0], ...
            'EdgeColor',       'none',        ...
            'FaceLighting',    'gouraud',     ...
            'AmbientStrength', 0.15);
view([-135 35]);
%NOTE: When rendering/loading in SolidWorks, use "point cloud" option 
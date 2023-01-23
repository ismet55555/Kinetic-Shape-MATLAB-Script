clc         %Clear the command window
close all   %Close all open figures/matlab windows
clear all   %Clear all data in matlab workspace

syms theta1  %Define Angle theta symbolic variable
syms LinVel1 %Define Linear velocity symbolic variable

%==============================================================================
%                                    INPUT 
%==============================================================================
%Angle Range Around Kinetic Shape                    %Unit Radians
Resolution  = 0.01;     %Discrete Devisions
Begin       = 0;        %Starting Angle
End         = 2*pi;     %Ending Angle
theta = [Begin : Resolution : End]';

%Initial Shape Definition Radius  
Ri = 2.75*0.0254;                                    %Unit Length

%Radial Input Force 
%     Note:  Has to be as a function of "theta1"    
Fr_Input = 5*cos(theta1) + 8;                       %Unit Force

%Applied Vertical Input Force (Applied Weight) 
%     Note:  Has to be as a function of "theta1"
Fv_Input = 25;                                      %Unit Force

%Intitial Angular Position and Angular Velocity (For Kinematics)
Initial_AngPos = 6.0;                                %Unit Radians
Initial_AngVel = 0;                                  %Unit Radians/Second

%Kinetic Shape Material (For Kinematics)
Density   = 632.68;                                  %Unit Mass/Length^3
Thickness = 0.0159;                                  %Unit Length

%Axle Mass or Dispensed Platform Mass (For Kinematics)
%  Note:  Enter '0' if none present
Mass2 = 5.5;                                         %Unit Mass

%Total Simulation Time and Time Resolution (For Kinematics)
TotalTime = 0.850;                                   %Unit Seconds
dt        = 0.0025;                                  %Unit Seconds

%After the Kinetic Shape is derieved, it is possible
%to apply another vertical applied force onto the shape
%and observe how it move. (For Kinematics)
%('theta1' as Variable) 
%    Note:  Fv_Input_App = Fv_Input to keep original applied force
Fv_Input_App = 14.5*4.45;                            %Unit Force

%Linear Rolling/Dispenser friction Coefficient Model (For Kinematics)
%  This is the resistance of linear motion 
%  between ground and axle
%('LinVel1' as Variable)
Mu_k = 0.145;                     %Kinetic Friction Coef   - Unitless
Mu_s = 0.450;                     %Static Friction Coef    - Unitless
C    = 8.0;                       %Friction Decay Constant - Unitless
Mu_L_fric = sign(LinVel1)*(Mu_k + (Mu_s - Mu_k)*exp(-abs(LinVel1)*C))...
           + Mu_k*LinVel1;                       %Unitless

%Over-ground rolling option
%Set Ground = 1 for when shape rolls over ground contact
%Set Ground = 0 for when shape rolls around its origin
Ground = 0;
%==============================================================================
%==============================================================================


fprintf('\n================================================\n')
fprintf('Initial Orientation Angle: %4.2f Radians\n',         Initial_AngPos)
fprintf('                           %4.2f Degrees\n',         Initial_AngPos*(180/pi))
fprintf('Initial Angular Velocity:  %4.2f Radians/Seconds\n', Initial_AngVel)
fprintf('                           %4.2f Degrees/Seconds\n', Initial_AngPos*(180/pi))
fprintf('================================================\n')

%Check if shape is contionous all around
%Output-print if continous or not
Shape_Cont = eval(int(Fr_Input, theta1, 0, 2*pi));
if Shape_Cont < 40
   fprintf('    ** Kinetic Shape is Continuous all Around. **\n')
else
  fprintf('     ** Kinetic Shape is NOT Continuous. **\n')
end
fprintf('================================================\n')

%Output-print Kinetic Shape Definition Forces
fprintf('Given/Desired Radial Force Input:\n') 
Fr_Input
fprintf('================================================\n')
fprintf('Applied Vertical Force Input:\n') 
Fv_Input
fprintf('================================================\n\n')



%==============================================================================
%               Applying Two-Dimensional Kinetic Shape Equation 
%==============================================================================
Rad = zeros(length(theta), 1);      %Creating Blank Vector
Rad(1) = Ri;                        %Initial Kinetic Shape Radius

Fr     = zeros(length(theta),1);    %Creating Blank Vectors
Fv     = zeros(length(theta),1);
Fv_app = zeros(length(theta),1);

F_Input = Fr_Input / Fv_Input;     %Devide Radial by Vertical Force
Fi      = int(F_Input, theta1);    %Inegrate Radial Force wrt theta1

ds   = 0;     
Area = 0;       
J    = 0;
for i = 2:length(theta) 
	%Subbing values into integral equation  
    F_int_1 = subs(Fi, [theta1],[theta(i)]);
    F_int_2 = subs(Fi, [theta1],[theta(i-1)]);
    
    %Derive Shape Radius, Rad(theta)
    Rad(i) = exp(  (F_int_1 - F_int_2)  ) * Rad(i-1);
    
    %Evaluating Specified Functions at Angle Step
    Fr(i)     = subs(Fr_Input,    [theta1],[theta(i)]);%Radial Force
    Fv(i)     = subs(Fv_Input,    [theta1],[theta(i)]);%Vertical Force
    Fv_app(i) = subs(Fv_Input_App,[theta1],[theta(i)]);%Vertical Force (Post)
    
    %Finding other Parameters
    dtheta = theta(i)-theta(i-1);                  %Step Angle (Radians)
    Area   = Area + (  (Rad(i)*dtheta)/2)*(Rad(i));%Shape Area (Length^2)
    ds     = ds + Rad(i) * dtheta;                 %Shape Arc Length (Length)
    J      = J + ( (Rad(i)^4)/4 ) * dtheta;        %Polar Moment of Inertia 
                                                   %      (Length^4)
end
Mass = Area*Density*Thickness;                     %Shape Mass

%==============================================================================
%                           Backward Check/Verification 
%==============================================================================
dr_dt_check     = zeros(length(Fr),1);    %Derivative of Shape Vector
Fr_check        = zeros(length(Fr),1);    %Horizontal Force Vector
%Comparing Given and Calculated Horizontal Reaction Force
for i = 2:length(Fr)
    %Derrivative of shape
    dr_dt_check(i) = (Rad(i)-Rad(i-1)) / (theta(i)-theta(i-1));
    
    %Polar Tangential Angle
    psi_check = atan( Rad(i) / dr_dt_check(i));
    
    %Horizontal Force From Derived Shape
    Fr_check(i) = Fv(i) * (cos(psi_check) / sin(psi_check));
end


%Command Window Parameter Output
fprintf('Minimum Shape Radius:\t %2.3f\t\t(Length)\n', Ri)
fprintf('Maximum Shape Radius:\t %2.3f\t\t(Length)\n', Rad(end))
fprintf('Shape Arc Length:\t\t %2.3f\t\t(Length)\n', ds)
fprintf('Shape Area:\t\t\t\t %2.3f\t\t(Length^2)\n', Area)
fprintf('Polar Moment of Inertia: %2.3f\t(Length^4)\n', J)
fprintf('Shape Mass:\t\t\t\t %2.3f\t\t(Mass)\n', Mass)
fprintf('=========================================\n\n') 


%==============================================================================
%             Mapping Angular Acceleration Around the Kinetic Shape    
%==============================================================================
dtheta     = mean( diff(theta) );      %Angle change/step
AngAcc     = zeros(length(theta), 1);  %Blank vector variables
J_ground   = zeros(length(theta), 1);
psi        = zeros(length(theta), 1);

for i = 2:length(Fr) 
    %Polar Mass Moment of Inertia for Shape (Mass*Length^2)
    J_mass    = Density * Thickness * J;        %Polar mass moment of inertia
    
    %Polar Mass Moment of Inertia For Axle or Dispenser (Mass*Length^2)
    J_Lin = Mass2*Rad(i)^2 * sin(psi(i));

    %Apply Parallel Axis Theorem for ground contact rolling
    J_ground(i) = J_mass + J_Lin + Ground*Mass*Rad(i)^2;  %Parallel Axis Thrm
    
    %Shape Radius Function Derrivative (Length/Radian)
    drdth = (Rad(i)-Rad(i-1)) / (theta(i)-theta(i-1)); 
    
    %Polar Tangential Angle (Radian)
    psi(i) = atan2( Rad(i), drdth );    
    
    %Angular Acceleration Value around Perimeter (Radian/Time^2)
    AngAcc(i) = Fv(i) * Rad(i)*cos( psi(i) ) / J_ground(i); 
end




%==============================================================================
%                     Simulating Kinetics Shape Movement with Time   
%==============================================================================
%Expanding Vectors
Rad2 = Rad;                             
for i = 1:4
   Fr       = [Fr, Fr];                 Fv       = [Fv, Fv];
   Fv_app   = [Fv_app, Fv_app];         Rad2     = [Rad2; Rad2];             
   psi      = [psi; psi];               drdth    = [drdth; drdth];
   J_ground = [J_ground; J_ground];
end
theta2 = [];                            %Expanding and Mirroring theta Vector
for i = 0:4;   theta2 = [theta2; theta+i*2*pi];   end
theta2 = [-theta2(end:-1:2); theta2];


time_index = linspace(1, TotalTime/dt, TotalTime/dt); %Creating Time Vector
time = zeros(length(time_index), 1);                  %Time Index Vector

%Blank vector variables
Mu_Rad_fric  = zeros(length(time_index), 1);    %Rotation Friction Model
Mu_Lin_fric  = zeros(length(time_index), 1);    %Linear Friction Model
AngPos       = zeros(length(time_index), 1);    %Angular Position
AngVel       = zeros(length(time_index), 1);    %Angular Velocity  
AngAcc2      = zeros(length(time_index), 1);    %Angular Acceleration  
LinAcc       = zeros(length(time_index), 1);    %Linear Position
LinVel       = zeros(length(time_index), 1);    %Linear Velocity
LinPos       = zeros(length(time_index), 1);    %Linear Acceleration 

%Initial Conditions (Initial Position, Velocity, Acceleration)
AngPos(1)    = Initial_AngPos;            
theta_tri    = delaunayn(theta2);                      %Triangulate 'theta'  
AngPos_idx   = dsearchn(theta2, theta_tri, AngPos(1));%Find initial angle index
                                                       
AngVel(1)   = Initial_AngVel;                          %Radians/Second
LinVel(1)   = AngVel(1)*Rad2(AngPos_idx)*sin(psi(1));  %Length/Second

Mu_Lin_fric(1) = -subs(Mu_L_fric, [LinVel1, theta1],...
                  [LinVel(1), theta2(AngPos_idx)]);         %Unitless
F_Lin_fric     = Fv_app(AngPos_idx)*Rad2(AngPos_idx)*sin(psi(AngPos_idx))...
                 *Mu_Lin_fric(1);                            %Force*Length

AngAcc2(1)  = -(  (Fv_app(AngPos_idx)*Rad2(AngPos_idx)*cos(psi(AngPos_idx))...
              - F_Lin_fric) / J_ground(AngPos_idx)  ) ; 

         
for t = time_index(2:end) 
    Mu_Lin_fric(t) = -subs(Mu_L_fric, [LinVel1, theta1],...
                      [LinVel(t-1), theta2(AngPos_idx)]);  %Unitless
    F_Lin_fric     = Fv_app(AngPos_idx)*Rad2(AngPos_idx)*sin(psi(AngPos_idx))...
                     *Mu_Lin_fric(t);                      %Force*Length
    
    %=========== Kinetic Shape Angular Kinetmatics ===========
    AngAcc2(t) = -((Fv_app(AngPos_idx)*Rad2(AngPos_idx)*cos(psi(AngPos_idx))...
                 - F_Lin_fric) / J_ground(AngPos_idx)); 
           
    if isnan( AngAcc2(t) ) == 1;    %Check if AngAcc2 is a number
        AngAcc2(t) = 0;             %If not, assign zero value
    end
    
    %Numerical Integration
    AngVel(t) = AngVel(t-1) + AngAcc2(t)*dt;   %Angular Velocity 
    AngPos(t) = AngPos(t-1) + AngVel(t) *dt;   %Angular Position 
        
    
    %============ Kinetic Shape Linear Kinematics =============
    LinAcc(t) = AngAcc2(t)* Rad2(AngPos_idx)*sin(psi(t));
    LinVel(t) = AngVel(t)* Rad2(AngPos_idx)*sin(psi(t));
    LinPos(t) = LinPos(t-1) + LinVel(t) *dt;   %Position(Numerical Integration)

    AngPos_idx = dsearchn(theta2, theta_tri, AngPos(t)); %Find angle index
    time(t) = time(t-1) + dt;                            %Increment Time Step
end


%==============================================================================
%              Plotting Kinetic Shape, Kinetics, and Kinematics  
%==============================================================================
%Initial Angle Index
theta_tri = delaunayn(theta);       %Triangulat 'theta'  
AngPos_idx = dsearchn(theta, theta_tri, Initial_AngPos); %Find angle index

figure(1);  %Plotting Wheel Shapes in Polar Coordinates
polar(theta, Rad, '.k')
hold on
h1 = polar(theta(AngPos_idx), Rad(AngPos_idx), 'ob');
set( findobj(h1, 'Type', 'line'), 'MarkerEdgeColor','b',...
     'MarkerFaceColor','b', 'MarkerSize',9);
title('Two Dimensional Kinetic Shape (Polar)','FontSize',15,...
      'FontName','Times New Roman')
xlabel('Theta (Degrees)','FontSize',12, 'FontName','Times New Roman')
ylabel('Radius (Inches)','FontSize',12, 'FontName','Times New Roman')
ylabh = get(gca,'ylabel');
set(ylabh,'Position',get(ylabh,'Position') - [.2 .2 0])
grid on

figure(2);  %Plotting Wheel Shapes in Cartasian Coordinates
hold on
plot(theta, Rad, '-k','Linewidth',2)
plot(theta(AngPos_idx), Rad(AngPos_idx),'ob','MarkerEdgeColor','b',...
     'MarkerFaceColor','b', 'MarkerSize',9);
yLimits = get(gca, 'YTick');
plot([theta(AngPos_idx), theta(AngPos_idx)], ...
     [min(yLimits), max(yLimits)],'--k', 'Linewidth', 2);
title('Two Dimensional Kinetic Shape (Cartesian)','FontSize',15,...
      'FontName','Times New Roman')
xlabel('Theta (Radians)','FontSize',12, 'FontName','Times New Roman')
ylabel('Radius (Inches)','FontSize',12, 'FontName','Times New Roman')
set(gca,'xLim',[0,2*pi]);
set(gca,'XTick',[0:pi/2:2*pi])
grid on


figure(3);  %Plotting Given and Calculated Horizontal Reaction Force
plot(theta, Fr, '--r','Linewidth',6)
hold on
plot(theta, Fv, '--b','Linewidth',6)
plot(theta, Fr_check, '-k','Linewidth',3)
plot([0 2*pi],[0 0],'k--','LineWidth',2)        %Zero Line
yLimits = get(gca, 'YTick');
plot([theta(AngPos_idx), theta(AngPos_idx)], [min(yLimits), max(yLimits)],...
      '--k', 'Linewidth', 2);  %Initial Line
title('Kinetic Shape Applied and Reaction Force','FontSize',15, 'FontName','Times New Roman')
xlabel('Theta (Radians)','FontSize',12, 'FontName','Times New Roman')
ylabel('Force (Newtons)','FontSize',12, 'FontName','Times New Roman')
set(gca,'xLim',[0,2*pi]);
set(gca,'XTick',[0:pi/2:2*pi])
grid on

figure(4)  %Kinetic Shape Angular Acceleration around Perimenter
plot(theta, AngAcc, '-.k', 'linewidth', 4)
hold on
yLimits = get(gca, 'YTick');
plot([theta(AngPos_idx), theta(AngPos_idx)], [min(yLimits), max(yLimits)], '--k', 'Linewidth', 2);
title({'Kinetic Shape Angular';'Acceleration at all Points'},'FontSize',15, 'FontName','Times New Roman')
set(gca,'xLim',[0,2*pi]);         set(gca,'XTick',[0:pi/2:2*pi])
xlabel('Theata(Radians)');        
ylabel('Angular Acceleration (Rad/s^2)','FontSize',12, 'FontName','Times New Roman')
grid on

figure(5) %KS Angular Position
plot(time, AngPos, '--g', 'linewidth', 3)
title('Angular Position','FontSize',15, 'FontName','Times New Roman')
xlabel('Time(Seconds)','FontSize',12, 'FontName','Times New Roman')
ylabel('Angle (Radians)','FontSize',12, 'FontName','Times New Roman')

figure(6) %KS Angular Velocity
plot(time, AngVel, '--b', 'linewidth', 3)
title('Angular Velocity','FontSize',15, 'FontName','Times New Roman')
xlabel('Time(Seconds)','FontSize',12, 'FontName','Times New Roman')
ylabel('Angular Velocity (Radians/Second)','FontSize',12, 'FontName','Times New Roman')

figure(7) %KS Angular Acceleration
plot(time, AngAcc2, '-.k', 'linewidth', 4)
title('Angular Acceleration','FontSize',15, 'FontName','Times New Roman')
xlabel('Time(Seconds)','FontSize',12, 'FontName','Times New Roman')
ylabel('Angular Acceleration (Radians/Second^2)','FontSize',12, 'FontName','Times New Roman')

figure(8) %KS Linear Position
plot(time, LinPos, '--g', 'linewidth', 3)
title({'Linear (Rolling) Position'},'FontSize',15, 'FontName','Times New Roman')
xlabel('Time(Seconds)','FontSize',12, 'FontName','Times New Roman')
ylabel('Position (Length)','FontSize',12, 'FontName','Times New Roman')

figure(9) %KS Linear Velocity
plot(time, LinVel, '--b', 'linewidth', 3)
title({'Linear (Rolling) Velocity'},'FontSize',15, 'FontName','Times New Roman')
xlabel('Time(Seconds)','FontSize',12, 'FontName','Times New Roman')
ylabel('Velocity (Length/Second)','FontSize',12, 'FontName','Times New Roman')

figure(10) %KS Linear Acceleration
plot(time, LinAcc, '-.k', 'linewidth', 4)
title({'Linear (Rolling) Acceleleration'},'FontSize',15, 'FontName','Times New Roman')
xlabel('Time(Seconds)','FontSize',12, 'FontName','Times New Roman')
ylabel('Acceleration (Length/Second^2)','FontSize',12, 'FontName','Times New Roman')

%==============================================================================
%                           Position and Format Plots  
%==============================================================================
%Positions of Plots on Screen
Pos  = zeros(10,4);
Pos(1,:)   = [0     0    440   348]; %KS in polar coordinates
Pos(2,:)   = [467   0    440   347]; %KS in Cartesian coordinates
Pos(3,:)   = [0     459  440   540]; %Applied/Reaction Forces
Pos(4,:)   = [470   459  440   540]; %Angular Acceleration around KS
Pos(5,:)   = [927   730  475   267]; %KS Angular Position
Pos(6,:)   = [927   379  475   258]; %KS Angular Velocity
Pos(7,:)   = [929   12   476   271]; %KS Angular Acceleration
Pos(8,:)   = [1420  730  480   267]; %KS Linear Position
Pos(9,:)   = [1420  377  480   263]; %KS Linear Velocity
Pos(10,:)  = [1421  16   481   268]; %KS Linear Acceleration

for i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    figure(i)
    set(gca, 'FontSize',10, 'FontName','Times New Roman')
    xlabel('Time(Seconds)')
    set(figure(i),'position', Pos(i,:)); 
    grid on
end

%==============================================================================
%                        Export 2D data points to Text File  
%==============================================================================
[X,Y] = pol2cart(theta, Rad);
Z = zeros(length(theta),1);
DATA = [X Y Z];
dlmwrite('2D_KS_XY.txt', DATA, 'precision', '%g', 'newline', 'pc');

%==============================================================================
%                       Export Dynamics Data to Text File  
%==============================================================================
DATA = [time, AngPos, AngVel, AngAcc2, LinPos, LinVel, LinAcc];
dlmwrite('2D_KS_KIN.txt', DATA, 'precision', '%g', 'newline', 'pc');

fprintf('Done.\n\n')

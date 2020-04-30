%% Biot-Savart law for a straight wire in x
% Plots a 3D numerical solution of the Biot-Savart integral
% for a straight wire in the x-dimension.
% C.L. James 2020
%----------------------------------------------------------%
clear;
close all;
clc;

% Physical constants
mu_0 = 1.25663706212e-6;
I    = 1;
c    = mu_0*I/(4*pi);
% Dimensions
L = 10;
n = L/2;
N = 1e3;
x = linspace(-n,n,N);
% y = linspace(-n,n,N);
% z = linspace(-n,n,N);

co = n+1; % coord-offset
b_y = zeros(L+1,L+1,L+1);
b_z = zeros(L+1,L+1,L+1);

for x_0 = -n:n
    for y_0 = -n:n
        for z_0 = -n:n
            
            for i = 1:length(x)-1
                
                x_seg_length = abs(x(i+1)-x(i));
                s_mid    = x(i)+x_seg_length/2;
                
                b1 =  (x_seg_length)/(((x_0-s_mid)^2 +y_0^2+z_0^2)^(3/2));                
                b_y(x_0+co,y_0+co,z_0+co) = b_y(x_0+co,y_0+co,z_0+co) ...
                    +z_0*b1;
                b_z(x_0+co,y_0+co,z_0+co) = b_z(x_0+co,y_0+co,z_0+co) ...
                    +y_0*b1;
                
            end
            
        end    
    end  
end

% Multiply loop output (b_y,b_z) by constants
B_y = -c*b_y;
B_z = c*b_z;
% Permute arrays so columns: x dimension, rows: y dimension
B_y = permute(B_y,[2,1,3]);
B_z = permute(B_z,[2,1,3]);
%Create 3d array of zeros for B_x
B_x = zeros(L+1,L+1,L+1);

% Create a meshgrid on which to plot data
xx = -n:n;
yy = -n:n;
zz = -n:n;
[X,Y,Z] = meshgrid(xx,yy,zz);

% Plot
quiver3(X,Y,Z,B_x, B_y, B_z)
xlabel('x');
ylabel('y');
zlabel('z');

% what is the field at (2,2,2)
ch1=B_y(2+co,2+co,2+co);
t_1 = atan((2+n)/sqrt(8));
t_2 = atan((2-n)/sqrt(8));
ch2 = 0.25*c*(sin(t_2)-sin(t_1));
% numerical = analytical iff ch1 = ch2


function LC_out = Generate_LCm_decoupled( n_x, n_y, A_size, lens_dist, adjacency )
%GENERATE_LCm_decoupled function initializes a set of light cones 
%associated with the sensor pixels behind one lenslet. The initial 
%angular span of the generated light cones is defined with respect to 
%the specified aperture.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   One light cone per pixel is generated. 
%   The centre of the image sensor has (0,0,0) coordinates.
%   The vertex of each light cone is located at the center of the pixel
%   area.

%   The input parameters are:
%   "n_x" is the number of pixels in x direction,
%   "n_y" is the number of pixels in y direction, 
%   "A_size" is the size of the aperture opening, 
%   "lens_dist" is the gap size between the sensor and the lensarray, 
%   "adjacency" is the distance between centres of two neighbouring pixels

%   The span of each light cone is determined by the size of the
%   aperture_opening

%   all distnaces are in um

%   n_x = 25 ;
%   n_y = 25
%   adjacency = 0.01;
%   acceptance_angle = 30/2/pi;


% clear all
% n_x = 2761;
% n_y = 2761;
% acceptance_angle = 1*2*pi/360;
% adjacency = 16.7;


%   ini_span = [teta_s, teta_f, phi_s, phi_f] in degrees, 
%   coming from the sensor's light-acceptance-angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LC_out = zeros(n_x * n_y,10);

for i = 1:n_x
    for j = 1:n_y
        LC_out((j-1) * n_x + i,1) = (-adjacency)*((n_x+1)/2) + (i * adjacency); % x_ini
        LC_out((j-1) * n_x + i,2) = (-adjacency)*((n_y+1)/2) + (j * adjacency); % y_ini
        LC_out((j-1) * n_x + i,4) = (-adjacency)*((n_x+1)/2) + (i * adjacency); % x
        LC_out((j-1) * n_x + i,5) = (-adjacency)*((n_y+1)/2) + (j * adjacency); % y
        LC_out((j-1) * n_x + i,7) = atan((-(A_size/2) - LC_out((j-1) * n_x + i,4)) / lens_dist); %span1
        LC_out((j-1) * n_x + i,8) = atan(((A_size/2) - LC_out((j-1) * n_x + i,4)) / lens_dist); %span2
        LC_out((j-1) * n_x + i,9) = atan((-(A_size/2) - LC_out((j-1) * n_x + i,5)) / lens_dist); %span3
        LC_out((j-1) * n_x + i,10) = atan(((A_size/2) - LC_out((j-1) * n_x + i,5)) / lens_dist); %span4
        
    end
end
LC_out(i,3) = 0; % z_ini
LC_out(i,6) = 0; % z
end


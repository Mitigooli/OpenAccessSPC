function LC_out = Generate_LCm( n_x, n_y,  acceptance_angle, adjacency )
%GENERATE_LCm function initializes a set of light cones
%associated with the sensor pixels.
%The angular span of the generated light cones is defined with respect to
%the given light-acceptance-angle of the sensor pixel as a physical
%property of the image sensor pixel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   One light cone per pixel is generated.
%   The centre of the image sensor has (0,0,0) coordinates.
%   The vertex of each light cone is located at the center of the pixel
%   area.

%   The input parameters are:
%   "n_x" is the number of pixels in x direction.
%   "n_y" is the number of pixels in y direction.
%   "acceptance_angle" is the light acceptance angle of each pixel in
%   degrees and with respect to the normal of the plane.
%   "adjacency" is the distance between centres of two neighbouring pixels

%   all distnaces are in um
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LC_out = zeros(n_x * n_y,10);

for i = 1:n_x
    for j = 1:n_y
        LC_out((j-1) * n_x + i,1) = (-adjacency)*((n_x+1)/2) + (i * adjacency); % x_ini
        LC_out((j-1) * n_x + i,2) = (-adjacency)*((n_y+1)/2) + (j * adjacency); % y_ini
        %       LC_out(i,3) = 0; % z_ini
        LC_out((j-1) * n_x + i,4) = (-adjacency)*((n_x+1)/2) + (i * adjacency); % x
        LC_out((j-1) * n_x + i,5) = (-adjacency)*((n_y+1)/2) + (j * adjacency); % y
        %       LC_out(i,6) = 0; % z
    end
end
LC_out(:,7) = - acceptance_angle;
LC_out(:,8) = acceptance_angle;
LC_out(:,9) = - acceptance_angle;
LC_out(:,10) = acceptance_angle;
end


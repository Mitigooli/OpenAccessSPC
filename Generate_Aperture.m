function A_out = Generate_Aperture( n_x, n_y, A_pitch, A_size, z )
%GENERATE_APERTURE will generate the aperture entity given the geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A_out is the generated aperture entity.
%   In case of a single aperture, A_out will have a single element.
%   In case of an aperture array, this function generates a rectangular 
%   grid of apertures, in which all the apertures have the same opening 
%   size "A_size" and uniform horizontal and vertical spacing "A_pitch" 
%   between their centres. 
%   "n_x" is the number of apertures in the aperture array in x direction.
%   "n_y" is the number of apertures in the aperture array in y direction.
%   "z" is the depth plane where the aperture entity is located.

%   Each aperture opening of A_out is a vector [x1, x2, y1, y2, z], where 
%   "x1" is the start point of the aperture opening in x direction,
%   "x2" is the end point of the aperture opening in x direction,
%   "y1" is the start point of the aperture opening in y direction,
%   "y2" is the end point of the aperture opening in y direction,
%   and "z" is the aperture plane (perpendicular to the optical axis).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xs = (-A_pitch) * round((n_x - 1) / 2) - (A_size/2):A_pitch:A_pitch * round((n_x - 1) / 2) - (A_size/2);
xf = ((-A_pitch) * round((n_x - 1) / 2) - (A_size/2):A_pitch:A_pitch * round((n_x - 1) / 2) - (A_size/2)) + A_size;
ys = (-A_pitch) * round((n_y - 1) / 2) - (A_size/2):A_pitch:A_pitch * round((n_y - 1) / 2) - (A_size/2);
yf = ((-A_pitch) * round((n_y - 1) / 2) - (A_size/2):A_pitch:A_pitch * round((n_y - 1) / 2) - (A_size/2)) + A_size;

count=1;
for i=1:n_x
    for j=0:n_y-1
        A_out(count).x1 = xs(1,i);
        A_out(count).x2 = xf(1,i);
        A_out(count).y1 = ys(1,j+1);
        A_out(count).y2 = yf(1,j+1);
        A_out(count).z = z;
        count = count+1;
    end
end


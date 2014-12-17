function L_out = Generate_Lens( n_x, n_y, L_pitch, z ,f)
%GENERATE_LENS will generate the lens entity given the geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   L_out is the generated lens entity.
%   In case of a single lens, L_out will have a single element.
%   In case of a lens array, this function generates a rectangular
%   grid of lenses, in which all the lenses have uniform horizontal and
%   vertical spacing "L_pitch" between their centres.
%   "n_x" is the number of lenses in the lens array in x direction.
%   "n_y" is the number of lenses in the lens array in y direction.
%   "z" is the depth plane where the lens entity is located.
%   "f" is the focal length of the lens.

%   Each lens in L_out is a structure with four arguments:
%   "x", "y", "z" are the coordinates of the optical center of the lens,
%   and "f" is the focal length of the lens.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [(-L_pitch) * round((n_x - 1) / 2) :L_pitch:L_pitch * round((n_x - 1) / 2)];
y = [(-L_pitch) * round((n_y - 1) / 2) :L_pitch:L_pitch * round((n_y - 1) / 2)];

count=1;
for i=1:n_x
    for j=0:n_y-1
        L_out(count).x = x(1,i);
        L_out(count).y = y(1,j+1);
        L_out(count).z = z;
        L_out(count).f = f;
        count = count+1;
    end
    
end


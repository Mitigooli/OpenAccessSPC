function LC_out = Lens_m( L_in, LC_in )
%LENS function simulates the effect of a physical lens on a set of
%input light cones.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Lens (which is considered to be infinitely wide and positioned normal
%   to the z axis) transforms a single or an array of input light cones
%   "LC_in" into a new (or a new set of) light cone(s) "LC_out".
%   The vertex and the angular span of the input light cone(s) are
%   transformed using the lens equation.

%   L_in is the input lens entity which can be a single or
%   an array of lenses on a single plane.
%   In case of a single lens, L_in has one element.
%   Each element in L_in is a structure with
%   "x", "y" , "z" as the coordinates of the optical center of the lens,
%   and "f" as the focal length of the lens

%   LC_in is the array if input light cone(s).
%   Each light cone is a 10*1 vector:
%   [x_ini, y_ini, z_ini, x, y, z, teta_s, teta_f, phi_s, phi_f] carrying
%   the coordinates of the initial light cone vertex [x_ini, y_ini, z_ini],
%   the coordinates of the current light cone vertex [x, y, z],
%   and the angular span of the light cone in x and y directions,
%   [teta_s, teta_f, phi_s, phi_f].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;

for i = 1:size(L_in,2)
    for j = 1:size(LC_in,1)
        
        if L_in(i).f > 0    % convex Lens
            
            B_temp = Base_m(LC_in(j,:),L_in(i).z);
            Lx = L_in(i).x ;
            Ly = L_in(i).y ;
            Lz = L_in(i).z ;
            Lf = L_in(i).f ;
            Bx1 = B_temp(1,1); % x1
            Bx2 = B_temp(1,2); % x2
            By1 = B_temp(1,3); % y1
            By2 = B_temp(1,4); % y2
            
            % modifying LCs according to the lens equation,
            count = count+1;
            
            % first we find z_out of the new light cone and then find the
            % intersection of the center ray and the z=z_out plane to find
            % the vertex of the new light cone.
            % The angular span is then calculated using this new vertex.
            
            LC_out(count,1) = LC_in(j,1); % x_ini
            LC_out(count,2) = LC_in(j,2); % y_ini
            LC_out(count,3) = LC_in(j,3); % z_ini
            LC_out(count,6) = Lz + ((Lz - LC_in(j,6)) * Lf) / (Lz - LC_in(j,6) - Lf ); % z
            LC_out(count,5) = LC_in(j,5)+((LC_in(j,6) - LC_out(count,6)) / (LC_in(j,6) - Lz)) * (Ly - LC_in(j,5)); % y
            LC_out(count,4) = LC_in(j,4)+((LC_in(j,6) - LC_out(count,6)) / (LC_in(j,6) - Lz)) * (Lx - LC_in(j,4)); % x
            LC_out(count,7) = atan((LC_out(count,4) - Bx1)/(LC_out(count,6) - Lz)); % span(1)
            LC_out(count,8) = atan((LC_out(count,4) - Bx2)/(LC_out(count,6) - Lz)); % span(2)
            LC_out(count,9) = atan((LC_out(count,5) - By1)/(LC_out(count,6) - Lz)); % span(3)
            LC_out(count,10) = atan((LC_out(count,5) - By2)/(LC_out(count,6) - Lz)); % span(4)
            
            % putting the new angular span in ascending order
            if LC_out(count,7) > LC_out(count,8)
                temp = LC_out(count,7);
                LC_out(count,7) = LC_out(count,8);
                LC_out(count,8) = temp;
            end
            if LC_out(count,9) > LC_out(count,10)
                temp = LC_out(count,9);
                LC_out(count,9) = LC_out(count,10);
                LC_out(count,10) = temp;
            end
        end
        
    end
end
if count == 0
    LC_out = 0;
end
end
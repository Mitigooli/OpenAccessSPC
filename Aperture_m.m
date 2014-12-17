function LC_out = Aperture_m( A_in, LC_in )
%APERTURE function simulates the effect of a physical aperture on a set of
%input light cones.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Aperture is applied on a single or an array of Light Cones (LC_in)
%   and generates a new (set of) Light cone(s) with trimmed span (LC_out).
%   For each input light cone, only the part of the light cone that can
%   pass the aperture opening will be present in the output light cone.

%   A_in is the input aperture entity which can be a single or
%   an array of apertures on a single plane.
%   In case of a single aperture, A_in has one element.
%   Each element in A_in is in the form of [x1, x2, y1, y2, z], where
%   x1 is the start point of the aperture opening in x direction,
%   x2 is the end point of the aperture opening in x direction,
%   y1 is the start point of the aperture opening in y direction,
%   y2 is the end point of the aperture opening in y direction,
%   and z is the aperture plane (perpendicular to the optical axis).

%   LC_in is the array of input light cone(s).
%   Each light cone is a 10*1 vector:
%   [x_ini, y_ini, z_ini, x, y, z, teta_s, teta_f, phi_s, phi_f] carrying
%   the coordinates of the initial light cone vertex [x_ini, y_ini, z_ini],
%   the coordinates of the current light cone vertex [x, y, z],
%   and the angular span of the light cone in x and y directions,
%   [teta_s, teta_f, phi_s, phi_f].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;

for i = 1:size(A_in,2)
    for j = 1:size(LC_in,1)
        % Finding the light-cone-base-area on the aperture plane and
        % the aperture opening
        B_temp = Base_m(LC_in(j,:),A_in(i).z);
        Ax1 = A_in(i).x1 ;
        Ax2 = A_in(i).x2 ;
        Ay1 = A_in(i).y1 ;
        Ay2 = A_in(i).y2 ;
        Bx1 = B_temp(1,1);
        Bx2 = B_temp(1,2);
        By1 = B_temp(1,3);
        By2 = B_temp(1,4);
        
        % modifying LCs if required
        C = [Ax1, Ay1, Ax2, Ay2];
        D = [Bx1, By1, Bx2, By2];
        if overlap( C, D ) > 0
            flag=1;
            sortx = sort([Ax1, Ax2, Bx1, Bx2]);
            sorty = sort([Ay1, Ay2, By1, By2]);
            count = count+1;
            LC_out(count,1) = LC_in(j,1); % x_ini
            LC_out(count,2) = LC_in(j,2); % y_ini
            LC_out(count,3) = LC_in(j,3); % z_ini
            LC_out(count,4) = LC_in(j,4); % x
            LC_out(count,5) = LC_in(j,5); % y
            LC_out(count,6) = LC_in(j,6); % z
            LC_out(count,7) = atan((sortx(1,2) - LC_in(j,4)) / (A_in(i).z - LC_in(j,6)));    % span(1)
            LC_out(count,8) = atan((sortx(1,3) - LC_in(j,4)) / (A_in(i).z - LC_in(j,6)));    % span(2)
            LC_out(count,9) = atan((sorty(1,2) - LC_in(j,5)) / (A_in(i).z - LC_in(j,6)));    % span(3)
            LC_out(count,10) = atan((sorty(1,3) - LC_in(j,5)) / (A_in(i).z - LC_in(j,6)));    % span(4)
            
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

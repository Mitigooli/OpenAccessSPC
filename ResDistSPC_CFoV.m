function res = ResDistSPC_CFoV (SPC, n_points, z_min, z_max)
%{
This function extracts the minimum resolvability distance from a known SPC 
at certain number of depth planes prependicular to the optical axis and in 
the Common-Field-of-View (CFoV).
"n_points" is the number of depth planes of interest.
"z_min" and "z_max" are maximum and minimum depth range.
--------------------------------------------------------------------------
Spatial resolution is defined as the inverse of the minimum resolvability
distance.
Resolvability distance is defined as the projected pixel size plus maximum
part size in x direction, and in the common field of view.
"z_plane" is the depth plane distance from the sensor plane
"pitch","n_l","g","n_x" are read along with SPC.
%}

if isempty(n_points)
    n_points = 1000;
end
if isempty(z_min)
    z_min = 250000;  %in um
end
if isempty(z_max)
    z_max = 3000000; %in um
end

load(SPC);

if z_min < n_l * g
    z_min = n_l * g;  % resolution planes > n_l * g , to be in the CFoV
end

ALCm = LC_tempL_m_final;
res = [];
z_plane = [];


for counter = 1/z_min:-((1/z_min-1/z_max)/n_points):1/z_max
    
    z = 1/counter;
    
    cfv = 2 * pitch / 2 /g * (z - n_l * g);
    Tb = Base_m(ALCm,z);
    
    for i = 1:size(Tb,1)
        Tb(i,13) = (Tb(i,1) + Tb(i,2)) / 2; % Finding the center of light cones
        Tb(i,14) = (Tb(i,3) + Tb(i,4)) / 2;
    end
    
    Tb = sortrows(Tb,[14 13 1]);
    temp = [NaN,NaN,NaN];
    
    % Finding "parts" (1D) in the common field of view
    for i = 1:size(Tb,1)-1
        if (Tb(i,14) == Tb(i+1,14))
            temp = [temp, Tb(i,1), Tb(i,2)];
        else
            if size(temp) ~= [1,3]
                break
            end
        end
    end
    
    
    for i = 1:size(temp,2)
        if (temp(i) > cfv/2) || (temp(i) < -cfv/2)
            temp(i) = NaN;
        end
    end
    
    temp = sort(temp);
    
    temp2 = temp(2:size(temp,2)) - temp(1:size(temp,2)-1);
    
    proj_pix_size = (z-g) / g * pitch / n_x;
    res = [res, proj_pix_size + max(abs(temp2))];
    z_plane = [z_plane,z];
end

figure;
hold on
plot(z_plane-g,res, 'r')

resFileName = ['ResDistSPC_CFoV_' num2str(n_points) 'points_' num2str(z_min) 'zmin_' num2str(z_max) 'zmax_' SPC];

save(resFileName,'z_plane','res','pitch','n_l','g','n_x','f','f_num','n_points','z_min','z_max');
end


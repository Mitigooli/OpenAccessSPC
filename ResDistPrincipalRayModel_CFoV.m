function res = ResDistPrincipalRayModel_CFoV (SPC, n_points, z_min, z_max)
%{
This function extracts the minimumresolvability distance using the 
principal-ray-model.
The resolvability distance is calculated for certain number of depth planes 
prependicular to the optical axis and only in the Common-Field-of-View, CFoV.
"n_points" is the number of depth planes of interest.
"z_min" and "z_max" are maximum and minimum depth range.
---------------------------------------------------------------------------
Spatial resolution is defined as the inverse of the minimum resolvability 
distance.
The resolvability distance is the maximun distance between the intersection 
of a principal ray and the depth plane of interest with the intersection of
the negiboring principal ray and the depth plane of interest.
--------------------------------------------------------------------------
Spatial resolution in the Principal Ray Model is obtained from the SPC model 
when each light cone is only the light ray at the center of the cone. 
In this case, the resolvability distance is obtained as: the maximum 
distance from the centre of the cone base on the depth plane of interest to 
its immediate neighbor in x direction, and in the common field of view.
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
        Tb(i,13) = (Tb(i,1) + Tb(i,2)) / 2;  % Principal ray position on the depth plane of interest
    end
    
    Tb = sortrows(Tb,[13 1]);
    temp = [0,0,0];
    
    
    % finding distances between the center point of each light-cone-base to its 
    % immediate neighbor in x direction, in the common field of view (CFoV)
    for i = 1:size(Tb,1)-1
        if (Tb(i,13) < cfv/2) && (Tb(i,13) > -cfv/2)
            Tb(i,15) = Tb(i+1,13) - Tb(i,13);
        else
            Tb(i,15) = 0;
        end;
    end
    
    res = [res,max(Tb(:,15))];
    z_plane = [z_plane,z];
end

figure;
hold on
plot(z_plane-g,2*res, 'r')

resFileName = ['ResDistPrincipalRayModel_CFoV_' num2str(n_points) 'points_' num2str(z_min) 'zmin_' num2str(z_max) 'zmax_' SPC];

save(resFileName,'z_plane','res','pitch','n_l','g','n_x','f','f_num','n_points','z_min','z_max')

end
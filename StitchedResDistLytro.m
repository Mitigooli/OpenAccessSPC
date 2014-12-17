% Stiched resolvability distance values for Lytro camera
% In each stitching stage, one lenslet is added.
% for each stitching stage, the SPC model of the subset of microlenses is
% generated and the resolvability distance is calculated for the CFoV of that
% subset of microlenses.
% The minimum distance from the microlens-plane is 2*g, 
% which is the plane where the CFoV of the first three leslets starts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lytro parameters:

f = 25;
f_num = 1.8;
pitch = 13.898614883422850808;
n_x = 10;
g = 25.1;
acceptance_angle = 40;

res = [];
z_plane = [];

for i = 328:-1:1
    
    n_l = i;
    
    aperture = f/f_num;                                                           % in um
    adjacency = pitch/n_x;                                                        % in um
    
    LCm = Generate_LCm(n_x, 1, acceptance_angle, adjacency );
    
    A = Generate_Aperture( 1, 1, pitch, aperture, g );
    
    microL = Generate_Lens( 1, 1, pitch, g, f);
    
    LC_tempA_m = Aperture_m( A, LCm);
    
    LC_tempL_m = Lens_m( microL, LC_tempA_m);  % Final LCs for a single lenslet
    
    LC_tempL_m_final = zeros(size(LC_tempL_m,1)*n_l , 10);
    
    for j = 1:n_l
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 1) = LC_tempL_m(:,1)+ pitch * (((n_l+1)/2)-j);       % x_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 2) = LC_tempL_m(:,2);                                % y_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 3) = LC_tempL_m(:,3);                                % z_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 4) = LC_tempL_m(:,4)+ pitch * (((n_l+1)/2)-j);       % x
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 5) = LC_tempL_m(:,5);                                % y
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 6) = LC_tempL_m(:,6);                                % z
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 7) = LC_tempL_m(:,7);                                % span1
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 8) = LC_tempL_m(:,8);                                % span2
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 9) = LC_tempL_m(:,9);                                % span3
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 10) = LC_tempL_m(:,10);                              % span4
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for zrange = -(i-1)*g  : g/20 : -(i-2)*g
        
        z = zrange;
        
        cfv = abs(pitch / g * (abs(z) - n_l * g));
        Tb = Base_m(LC_tempL_m_final,z);
        
        for j = 1:size(Tb,1)
            Tb(j,13) = (Tb(j,1) + Tb(j,2)) / 2;
            Tb(j,14) = (Tb(j,3) + Tb(j,4)) / 2;
        end
        
        Tb = sortrows(Tb,[14 13 1]);
        temp = [NaN,NaN,NaN];
        
        % finding start and end points in a horizontal line in the common field of view
        for j = 1:size(Tb,1)-1
            if (Tb(j,14) == Tb(j+1,14)) %&& (Tb(i,13) < cfv/2) && (Tb(i,13) > -cfv/2)
                temp = [temp, Tb(j,1), Tb(j,2)];
            else
                if size(temp) ~= [1,3]
                    break
                end
            end
        end
        
        
        for j = 1:size(temp,2)
            if (temp(j) > (cfv/2) + (pitch/2)) || (temp(j) < -((cfv/2) + (pitch/2)))
                temp(j) = NaN;
            end
        end
        
        temp = sort(temp);
        
        temp2 = temp(2:size(temp,2)) - temp(1:size(temp,2)-1);
        
        proj_pix_size = abs((z-g) / g * pitch / n_x);
        res = [res, proj_pix_size + max(abs(temp2))];
        z_plane = [z_plane,z];
    end
end

for i = 1:328
    
    n_l = i;
    
    aperture = f/f_num;                                                           % in um
    adjacency = pitch/n_x;                                                        % in um
    
    LCm = Generate_LCm(n_x, 1, acceptance_angle, adjacency );
    
    A = Generate_Aperture( 1, 1, pitch, aperture, g );
    
    microL = Generate_Lens( 1, 1, pitch, g, f);
    
    LC_tempA_m = Aperture_m( A, LC_tempA_m);
    
    LC_tempL_m = Lens_m( microL, LCm);  % Final LCs for a single lenslet
    
    LC_tempL_m_final = zeros(size(LC_tempL_m,1)*n_l , 10);
    
    for j = 1:n_l
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 1) = LC_tempL_m(:,1)+ pitch * (((n_l+1)/2)-j);       % x_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 2) = LC_tempL_m(:,2);                                % y_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 3) = LC_tempL_m(:,3);                                % z_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 4) = LC_tempL_m(:,4)+ pitch * (((n_l+1)/2)-j);       % x
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 5) = LC_tempL_m(:,5);                                % y
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 6) = LC_tempL_m(:,6);                                % z
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 7) = LC_tempL_m(:,7);                                % span1
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 8) = LC_tempL_m(:,8);                                % span2
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 9) = LC_tempL_m(:,9);                                % span3
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 10) = LC_tempL_m(:,10);                              % span4
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for zrange = i*g  : g/20 : (i+1)*g
        
        z = zrange;
        
        cfv = 2 * pitch / 2 /g * (z - n_l * g);
        Tb = Base_m(LC_tempL_m_final,z);
        
        for j = 1:size(Tb,1)
            Tb(j,13) = (Tb(j,1) + Tb(j,2)) / 2;
            Tb(j,14) = (Tb(j,3) + Tb(j,4)) / 2;
        end
        
        Tb = sortrows(Tb,[14 13 1]);
        temp = [NaN,NaN,NaN];
        
        % finding start and end points in a horizontal line in the common field of view
        for j = 1:size(Tb,1)-1
            if (Tb(j,14) == Tb(j+1,14)) %&& (Tb(i,13) < cfv/2) && (Tb(i,13) > -cfv/2)
                temp = [temp, Tb(j,1), Tb(j,2)];
            else
                if size(temp) ~= [1,3]
                    break
                end
            end
        end
        
        
        for j = 1:size(temp,2)
            if (temp(j) > (cfv/2) + (pitch/2)) || (temp(j) < -((cfv/2) + (pitch/2)))
                temp(j) = NaN;
            end
        end
        
        temp = sort(temp);
        
            
        temp2 = temp(2:size(temp,2)) - temp(1:size(temp,2)-1);
        
        proj_pix_size = (z-g) / g * pitch / n_x;
        res = [res, proj_pix_size + max(abs(temp2))];
        z_plane = [z_plane,z];
    end
end
i = 328;
for zrange = (i+1)*g + g/20 : g : 30000 ;
    
    z = zrange;
    
    cfv = 2 * pitch / 2 /g * (z - n_l * g);
    Tb = Base_m(LC_tempL_m_final,z);
    
    for j = 1:size(Tb,1)
        Tb(j,13) = (Tb(j,1) + Tb(j,2)) / 2;
        Tb(j,14) = (Tb(j,3) + Tb(j,4)) / 2;
    end
    
    Tb = sortrows(Tb,[14 13 1]);
    temp = [0,0,0];
    
    % finding "parts" (1D) in the common field of view
    for j = 1:size(Tb,1)-1
        if (Tb(j,14) == Tb(j+1,14))
            temp = [temp, Tb(j,1), Tb(j,2)];
        else
            if size(temp) ~= [1,3]
                break
            end
        end
    end
    
    
    for j = 1:size(temp,2)
        if (temp(j) > cfv/2) || (temp(j) < -cfv/2)
            temp(j) = 0;
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
plot(z_plane,res, 'r')

resFileName = 'StitchedResDistLytro';
save(resFileName,'z_plane','res');


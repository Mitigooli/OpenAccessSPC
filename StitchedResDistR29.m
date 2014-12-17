% Stiched minimum resolvability distance values for R29 camera
% In each stitching stage, three lenslets are added.
% for each stitching stage, the SPC model of the subset of microlenses is
% generated and the resolvability distance is calculated for the CFoV of that
% subset of microlenses.
% The minimum distance from the microlens-plane is 2*g,
% which is the plane where the CFoV of the first three leslets starts.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Reading the input regarding the camera geometry
%--------------------------------------------------------------------------

f1 = 1000*input('Focal length of the type I lenslet, f1(mm): ');              % in um
f2 = 1000*input('Focal length of the type I lenslet, f2(mm): ');              % in um
f3 = 1000*input('Focal length of the type I lenslet, f3(mm): ');              % in um
g = 1000*input('The gap between the sensor and the lenslet array, g(mm): ');  % in um

wf_num = input('working f-number of the lenslets, f#: ');
if isempty(wf_num)
    wf_num = 7;
end

pitch = 1000*input('Pitch size of the lenslet, p(mm): ');                     % in um
if isempty(pitch)
    pitch = 171;
end

n_l = input('Number of lenslets in a row: ');
if isempty(n_l)
    n_l = 210;
end

n_x = input('Number of pixels in a row behind one lenslet: ');
if isempty(n_x)
    n_x = 31;
end

acceptance_angle = 40;

res = [];
z_plane = [];

%Generating the SPC model for smaller number of lenslets
%and looking into the resolution in the CFoV in each case
for k = (n_l-3):-3:3
    
    n_l_temp = k;
    
    aperture = g/wf_num;                                                          % in um
    adjacency = pitch/n_x;                                                        % in um
    
    LCm = Generate_LCm(n_x, 1, acceptance_angle, adjacency );
    
    A = Generate_Aperture( 1, 1, pitch, aperture, g );
    
    microL1 = Generate_Lens( 1, 1, pitch, g, f1);
    microL2 = Generate_Lens( 1, 1, pitch, g, f2);
    microL3 = Generate_Lens( 1, 1, pitch, g, f3);
    
    LC_tempA_m = Aperture_m( A, LCm);
    
    LC_tempL_m1 = Lens_m( microL1, LC_tempA_m);  % Final LCs for a single lenslet (type I)
    LC_tempL_m2 = Lens_m( microL2, LC_tempA_m);  % Final LCs for a single lenslet (type II)
    LC_tempL_m3 = Lens_m( microL3, LC_tempA_m);  % Final LCs for a single lenslet (type III)
    
    LC_tempL_m_final = zeros(size(LC_tempL_m1,1)*n_l_temp , 10);
    
    for i = 1:3:n_l_temp
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 1) = LC_tempL_m1(:,1)+ pitch * (((n_l_temp+1)/2)-i);  % x_ini
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 2) = LC_tempL_m1(:,2);                                % y_ini
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 3) = LC_tempL_m1(:,3);                                % z_ini
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 4) = LC_tempL_m1(:,4)+ pitch * (((n_l_temp+1)/2)-i);  % x
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 5) = LC_tempL_m1(:,5);                                % y
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 6) = LC_tempL_m1(:,6);                                % z
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 7) = LC_tempL_m1(:,7);                                % span1
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 8) = LC_tempL_m1(:,8);                                % span2
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 9) = LC_tempL_m1(:,9);                                % span3
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 10) = LC_tempL_m1(:,10);                              % span4
        
        j = i+1;
        
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 1) = LC_tempL_m2(:,1)+ pitch * (((n_l_temp+1)/2)-j);  % x_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 2) = LC_tempL_m2(:,2);                                % y_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 3) = LC_tempL_m2(:,3);                                % z_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 4) = LC_tempL_m2(:,4)+ pitch * (((n_l_temp+1)/2)-j);  % x
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 5) = LC_tempL_m2(:,5);                                % y
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 6) = LC_tempL_m2(:,6);                                % z
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 7) = LC_tempL_m2(:,7);                                % span1
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 8) = LC_tempL_m2(:,8);                                % span2
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 9) = LC_tempL_m2(:,9);                                % span3
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 10) = LC_tempL_m2(:,10);                              % span4
        
        j = i+2;
        
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 1) = LC_tempL_m3(:,1)+ pitch * (((n_l_temp+1)/2)-j);  % x_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 2) = LC_tempL_m3(:,2);                                % y_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 3) = LC_tempL_m3(:,3);                                % z_ini
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 4) = LC_tempL_m3(:,4)+ pitch * (((n_l_temp+1)/2)-j);  % x
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 5) = LC_tempL_m3(:,5);                                % y
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 6) = LC_tempL_m3(:,6);                                % z
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 7) = LC_tempL_m3(:,7);                                % span1
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 8) = LC_tempL_m3(:,8);                                % span2
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 9) = LC_tempL_m3(:,9);                                % span3
        LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 10) = LC_tempL_m3(:,10);                              % span4
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for zrange = -(k+2)*g  : g/20 : -(k-1)*g % givig the depth planes of interest
        
        z = zrange;
        
        cfv = abs(pitch / g * (abs(z) - (n_l_temp-1) * g));
        Tb = Base_m(LC_tempL_m_final,z);
        
        for j = 1:size(Tb,1)
            Tb(j,13) = (Tb(j,1) + Tb(j,2)) / 2;
            Tb(j,14) = (Tb(j,3) + Tb(j,4)) / 2;
        end
        
        Tb = sortrows(Tb,[14 13 1]);
        temp = [NaN,NaN,NaN];
        
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
            if (temp(j) > (cfv/2)) || (temp(j) < -(cfv/2))
                temp(j) = NaN;
            end
        end
        
        temp = sort(temp);
        
        
        temp2 = temp(2:size(temp,2)) - temp(1:size(temp,2)-1);
        
        proj_pix_size = (abs(z)+g) / g * pitch / n_x;
        res = [res, proj_pix_size + max(abs(temp2))];
        z_plane = [z_plane,z];
    end
end

% uncomment if you want a plot
% figure;
% hold on
% plot(z_plane,res, 'b')

resFileName = 'StitchedResDistR29';
save(resFileName,'z_plane','res');


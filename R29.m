% Implementing the SPC model for a Single ROW of the sensor pixels in
% an R29 RAYTRIX CAMERA.
% Saving the corresponding SPC model in a .mat file

% Assumptions:
% Decoupled lenslets,
% Excluding the main objective lens,
% The lenslet's aperture opening size "aperture" and
% the lenslet's pitch "pitch" can be different.

% Default setup:
% Removed XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% The working f-number is considered, which is equal to
% g/aperture_opening_diameter and so is constant for all three types
% of microlenses.

%--------------------------------------------------------------------------
% Reading the input regarding the camera geometry
%--------------------------------------------------------------------------

f1 = 1000*input('Focal length of the type I lenslet, f1(mm): ');              % in um
f2 = 1000*input('Focal length of the type I lenslet, f2(mm): ');              % in um
f3 = 1000*input('Focal length of the type I lenslet, f3(mm): ');              % in um

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

g = 1000*input('The gap between the sensor and the lenslet array, g(mm): ');  % in um


aperture = g/wf_num;                                                          % in um
adjacency = pitch/n_x;                                                        % in um

if aperture > pitch
    
    message = ['g/f# should be smaller than the pitch of the lenslet'];
    
    pitch = 1000*input('Pitch size of the lenslet, p(mm): ');                 % in um
    if isempty(pitch)
        pitch = 171;
    end
    
    g = 1000*input('The gap between the sensor and the lenslet array, g(mm): ');  % in um
    
    
    wf_num = input('working f-number of the lenslets, f#: ');
    if isempty(wf_num)
        wf_num = 7;
    end
    
    aperture = g/wf_num;                                                      % in um
    adjacency = pitch/n_x;                                                    % in um
    
end
acceptance_angle = 40; %assumed pixel's light acceptance angle
%-----------------------------------------------------------------------------
% Initializing the Light Cones on the sensor array
%-----------------------------------------------------------------------------
% Note: here we implement a single lenslet setup!
% Later we build the multi lenslet setup as a union of single lenslet setups
% LC_out = Generate_LCm( n_x, n_y, acceptance_angle, adjacency )

LCm = Generate_LCm(n_x, 1, acceptance_angle, adjacency );

%-----------------------------------------------------------------------------
% initializing the single Aperture regarding one microlens
%-----------------------------------------------------------------------------
% A_out = Generate_Aperture( n_x, n_y, A_pitch, A_size, z )
% 1x1 test: A = Generate_Aperture( 1, 1, 4200, 1000*f/f#, 12510 );

A = Generate_Aperture( 1, 1, pitch, aperture, g );

%-----------------------------------------------------------------------------
% initializing the single Lens regarding one microlens (for each microles type)
%-----------------------------------------------------------------------------
% L_out = Generate_Lens( n_x, n_y, L_pitch, z ,f)

microL1 = Generate_Lens( 1, 1, pitch, g, f1);
microL2 = Generate_Lens( 1, 1, pitch, g, f2);
microL3 = Generate_Lens( 1, 1, pitch, g, f3);

%-----------------------------------------------------------------------------
% applying one aperture on all the light cones
%-----------------------------------------------------------------------------
% LC_out_m = Aperture_m(A_in, LCm_in )

LC_tempA_m = Aperture_m( A, LCm);

%-----------------------------------------------------------------------------
% applying the lens effect
%-----------------------------------------------------------------------------
% LC_out_m = Lens_m( L_in, LCm_in )

LC_tempL_m1 = Lens_m( microL1, LC_tempA_m);  % Final LCs for a single lenslet (type I)
LC_tempL_m2 = Lens_m( microL2, LC_tempA_m);  % Final LCs for a single lenslet (type II)
LC_tempL_m3 = Lens_m( microL3, LC_tempA_m);  % Final LCs for a single lenslet (type III)

%-----------------------------------------------------------------------------
% the SPC corresponding to an array of n_lx1 lenslets is generated
%-----------------------------------------------------------------------------

LC_tempL_m_final = zeros(size(LC_tempL_m1,1)*n_l , 10);

for i = 1:3:n_l
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 1) = LC_tempL_m1(:,1)+ pitch * (((n_l+1)/2)-i);       % x_ini
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 2) = LC_tempL_m1(:,2);                                % y_ini
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 3) = LC_tempL_m1(:,3);                                % z_ini
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 4) = LC_tempL_m1(:,4)+ pitch * (((n_l+1)/2)-i);       % x
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 5) = LC_tempL_m1(:,5);                                % y
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 6) = LC_tempL_m1(:,6);                                % z
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 7) = LC_tempL_m1(:,7);                                % span1
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 8) = LC_tempL_m1(:,8);                                % span2
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 9) = LC_tempL_m1(:,9);                                % span3
    LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 10) = LC_tempL_m1(:,10);                              % span4
    
    j = i+1;
    
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 1) = LC_tempL_m2(:,1)+ pitch * (((n_l+1)/2)-j);       % x_ini
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 2) = LC_tempL_m2(:,2);                                % y_ini
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 3) = LC_tempL_m2(:,3);                                % z_ini
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 4) = LC_tempL_m2(:,4)+ pitch * (((n_l+1)/2)-j);       % x
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 5) = LC_tempL_m2(:,5);                                % y
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 6) = LC_tempL_m2(:,6);                                % z
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 7) = LC_tempL_m2(:,7);                                % span1
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 8) = LC_tempL_m2(:,8);                                % span2
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 9) = LC_tempL_m2(:,9);                                % span3
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 10) = LC_tempL_m2(:,10);                              % span4
    
    j = i+2;
    
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 1) = LC_tempL_m3(:,1)+ pitch * (((n_l+1)/2)-j);       % x_ini
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 2) = LC_tempL_m3(:,2);                                % y_ini
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 3) = LC_tempL_m3(:,3);                                % z_ini
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 4) = LC_tempL_m3(:,4)+ pitch * (((n_l+1)/2)-j);       % x
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 5) = LC_tempL_m3(:,5);                                % y
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 6) = LC_tempL_m3(:,6);                                % z
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 7) = LC_tempL_m3(:,7);                                % span1
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 8) = LC_tempL_m3(:,8);                                % span2
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 9) = LC_tempL_m3(:,9);                                % span3
    LC_tempL_m_final(((j-1) * n_x) + 1 : j * n_x , 10) = LC_tempL_m3(:,10);                              % span4
    
end

%-----------------------------------------------------------------------------
% Saving the resulting sets of light cones geometry in the output files (Single lens and Multi lens formats)
%-----------------------------------------------------------------------------

% SingleLensletImplementation = ['1x1_lenslet_output_decoupled_' num2str(n_x) 'pixel_p=' num2str(pitch/1000) 'mm_g=' num2str(g/1000) 'mm_f=' num2str(f/1000) 'mm_f#=' num2str(f_num) '.mat'];
% save (SingleLensletImplementation,'LC_tempL_m','pitch','n_l','g','n_x');

MultiLensletImplementation = ['R29_' num2str(n_l) 'x1_lenslet_output_decoupled_' num2str(n_x) 'pixel_p=' num2str(pitch/1000) 'mm_g=' num2str(g/1000) 'mm_f1=' num2str(f1/1000) 'mm_f2=' num2str(f2/1000) 'mm_f3=' num2str(f3/1000) 'mm_wf#=' num2str(wf_num) '.mat'];
save (MultiLensletImplementation,'LC_tempL_m_final','pitch','n_l','g','n_x','f1','f2','f3','wf_num');



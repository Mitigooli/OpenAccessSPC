% Implementing the SPC model for a SINGLE ROW of the sensor pixels in 
% a LYTRO CAMERA.
% Saving the corresponding SPC model in a .mat file 

% Assumptions:
% decoupled lenslets,
% excluding the main objective lens,
% lenslet's aperture size "aperture" and the lenslet's pitch "pitch" can be different

% Default setup 1: 
% (328X1 lenslets, 3280 px, pitch=1.3898614883422850808e-02 mm, g=2.51-02 mm, f=2.5e-02 mm, f#=1.8)

%-----------------------------------------------------------------------------
% Reading the input regarding the camera geometry
%-----------------------------------------------------------------------------

f = 1000*input('Focal length of the lenslets, f(mm): ');                      % in um
if isempty(f)
    f = 25;
end

f_num = input('f-number of the lenslets, f#: ');
if isempty(f_num)
    f_num = 1.8;
end

pitch = 1000*input('Pitch size of the lenslet, p(mm): ');                     % in um
if isempty(pitch)
    pitch = 13.898614883422850808;
end

n_l = input('Number of lenslets in a row (an odd number): ');
if isempty(n_l)
    n_l = 327;
end

n_x = input('Number of pixels in a row behind one lenslet: ');
if isempty(n_x)
    n_x = 10;
end

g = 1000*input('The gap between the sensor and the lenslet array, g(mm): ');  % in um
if isempty(g)
    g = 26.315;
end

aperture = f/f_num;                                                           % in um
adjacency = pitch/n_x;                                                        % in um

if aperture > pitch
    
    message = ['f/f# should be smaller than the pitch of the lenslet'];
    
    pitch = 1000*input('Pitch size of the lenslet, p(mm): ');                 % in um
    if isempty(pitch)
        pitch = 13.898614883422850808;
    end
    
    f = 1000*input('Focal length of the lenslets, f(mm): ');                  % in um
    if isempty(f)
        f = 25;
    end
    
    f_num = input('f-number of the lenslets, f#: ');
    if isempty(f_num)
        f_num = 1.8;
    end
    
    aperture = f/f_num;                                                       % in um
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
% initializing the single Lens regarding one microlens
%-----------------------------------------------------------------------------
% L_out = Generate_Lens( n_x, n_y, L_pitch, z ,f)

microL = Generate_Lens( 1, 1, pitch, g, f);

%-----------------------------------------------------------------------------
% applying one aperture on all the light cones
%-----------------------------------------------------------------------------
% LC_out_m = Aperture_m(A_in, LCm_in )

LC_tempA_m = Aperture_m( A, LCm);

%-----------------------------------------------------------------------------
% applying the lens effect 
%-----------------------------------------------------------------------------
% LC_out_m = Lens_m( L_in, LCm_in )

LC_tempL_m = Lens_m( microL, LC_tempA_m);  % Final LCs for a single lenslet

%-----------------------------------------------------------------------------
% the data set corresponding to an array of n_lx1 lenslets is generated
%-----------------------------------------------------------------------------

LC_tempL_m_final = zeros(size(LC_tempL_m,1)*n_l , 10);

for i = 1:n_l
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 1) = LC_tempL_m(:,1)+ pitch * (((n_l+1)/2)-i);       % x_ini
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 2) = LC_tempL_m(:,2);                                % y_ini
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 3) = LC_tempL_m(:,3);                                % z_ini
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 4) = LC_tempL_m(:,4)+ pitch * (((n_l+1)/2)-i);       % x
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 5) = LC_tempL_m(:,5);                                % y
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 6) = LC_tempL_m(:,6);                                % z
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 7) = LC_tempL_m(:,7);                                % span1
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 8) = LC_tempL_m(:,8);                                % span2
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 9) = LC_tempL_m(:,9);                                % span3
        LC_tempL_m_final(((i-1) * n_x) + 1 : i * n_x , 10) = LC_tempL_m(:,10);                              % span4
end

%-----------------------------------------------------------------------------
% Saving the resulting sets of light cones geometry in the output files (Single lens and Multi lens formats)
%-----------------------------------------------------------------------------

% SingleLensletImplementation = ['1x1_lenslet_output_decoupled_' num2str(n_x) 'pixel_p=' num2str(pitch/1000) 'mm_g=' num2str(g/1000) 'mm_f=' num2str(f/1000) 'mm_f#=' num2str(f_num) '.mat'];
% save (SingleLensletImplementation,'LC_tempL_m','pitch','n_l','g','n_x');

MultiLensletImplementation = ['Lytro_' num2str(n_l) 'x1_lenslet_output_decoupled_' num2str(n_x) 'pixel_p=' num2str(pitch/1000) 'mm_g=' num2str(g/1000) 'mm_f=' num2str(f/1000) 'mm_f#=' num2str(f_num) '.mat'];
save (MultiLensletImplementation,'LC_tempL_m_final','pitch','n_l','g','n_x','f','f_num');


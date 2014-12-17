function Bm = Base_m( LCm, z0 )
%BASE provides the base of the Light Cone(s) when projected on the plane z=z0
%prependicular to the optical axis
%   Detailed explanation goes here


    Bm(:,1) = LCm(:,4)+(z0 - LCm(:,6)).*tan(LCm(:,7));  % x1
    Bm(:,2) = LCm(:,4)+(z0 - LCm(:,6)).*tan(LCm(:,8));  % x2
    Bm(:,3) = LCm(:,5)+(z0 - LCm(:,6)).*tan(LCm(:,9));  % y1
    Bm(:,4) = LCm(:,5)+(z0 - LCm(:,6)).*tan(LCm(:,10)); % y2
    Bm(:,5) = z0;
    Bm(:,6) = LCm(:,7); % span1
    Bm(:,7) = LCm(:,8); % span2
    Bm(:,8) = LCm(:,9); % span3
    Bm(:,9) = LCm(:,10);% span4
    Bm(:,10) = LCm(:,1);% initip1
    Bm(:,11) = LCm(:,2);% initip2
    Bm(:,12) = LCm(:,3);% initip3
%     Bm(i,1) = LCm(i,4)+(z0 - LCm(i,6))*tan(LCm(i,7));
%     Bm(i,2) = LCm(i,4)+(z0 - LCm(i,6))*tan(LCm(i,8));
%     Bm(i,3) = LCm(i,5)+(z0 - LCm(i,6))*tan(LCm(i,9));
%     Bm(i,4) = LCm(i,5)+(z0 - LCm(i,6))*tan(LCm(i,10));
%     Bm(i,5) = z0;
%     Bm(i,6) = LCm(i,7);
%     Bm(i,7) = LCm(i,8);
%     Bm(i,8) = LCm(i,9);
%     Bm(i,9) = LCm(i,10);
%     Bm(i,10) = LCm(i,1);
%     Bm(i,11) = LCm(i,2);
%     Bm(i,12) = LCm(i,3);
    
for i=1:size(LCm,1)    
    if Bm(i,1) > Bm(i,2)
        temp = Bm(i,1);
        Bm(i,1) = Bm(i,2);
        Bm(i,2) = temp;
    end
    if Bm(i,3) > Bm(i,4)
        temp = Bm(i,3);
        Bm(i,3) = Bm(i,4);
        Bm(i,4) = temp;
    end
end


% Plot_Base;

end


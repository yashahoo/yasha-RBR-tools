function [Est,Nth] = pol2cartCOMPASS(Theta, Speed)
% pol2cartCOMPASS takes in velocities in compass coordinates where Theta is
% the angle in compass coordinates (ie direction goes from North=0
% clockwise so East =90). returns velocities in East (x) and North (y)
% componants
%
%   [Est,Nth] = pol2cartCOMPASS(Theta, Speed)
%
%   Inputs:     Theta = angle of maximum variance, geographical notation
%               (north = 0, east = 90)
%               Speed = speed in user defined units of measurements
%
%   Output:     Est,Nth = velocities in the East and North direction
%               respectively in the same units as speed
%
%
% Author(s): SDG
% Created: 05-Dec-2007
% $Revision: 1 $  $Date: 2007/11/27 14:06:02 $
Theta=mod(Theta,360);
Est=zeros(size(Theta,1),size(Theta,2));
Nth=zeros(size(Theta,1),size(Theta,2));
for i=1:size(Theta,1)
    for j=1:size(Theta,2)
        if   Theta(i,j)==0
            Est(i,j)=0;
            Nth(i,j)=Speed(i,j);
        elseif Theta(i,j)<90
            Est(i,j)=Speed(i,j)*sind(Theta(i,j));
            Nth(i,j)=Speed(i,j)*cosd(Theta(i,j));
        elseif Theta(i,j)==90
            Est(i,j)=Speed(i,j);
            Nth(i,j)=0;
        elseif   Theta(i,j)<180
            Est(i,j)=Speed(i,j)*cosd(Theta(i,j)-90);
            Nth(i,j)=-Speed(i,j)*sind(Theta(i,j)-90);
        elseif   Theta(i,j)==180
            Est(i,j)=0;
            Nth(i,j)=-Speed(i,j);
        elseif   Theta(i,j)<270
            Est(i,j)=-Speed(i,j)*sind(Theta(i,j)-180);
            Nth(i,j)=-Speed(i,j)*cosd(Theta(i,j)-180);
        elseif   Theta(i,j)==270
            Est(i,j)=-Speed(i,j);
            Nth(i,j)=0;
        elseif   Theta(i,j)<360
            Est(i,j)=-Speed(i,j)*cosd(Theta(i,j)-270);
            Nth(i,j)=Speed(i,j)*sind(Theta(i,j)-270);
        elseif   Theta(i,j)==360
            Est(i,j)=0;
            Nth(i,j)=Speed(i,j);
        end
    end
end

return;
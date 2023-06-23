function [Theta Speed] = cart2polCOMPASS(Est,Nth)
% cart2polCOMPASS takes in velocities in the East (x) and North (y) direction and
% returns the same information in Direction Speed. Direction is returned in
% compass coordinates and in degrees
%
%   [Theta Speed] = cart2polCOMPASS(Est,Nth)
%
%   Inputs:     Est,Nth = velocities in the East and North direction
%               respectively Matricees can be 1 or 2 D but must be the same size   
%               
%
%   Output:     Theta = angle of maximum variance, geographical notation
%               (north = 0, east = 90)
%               Speed = sqrt(Est.^2+Nth.^2);
%
%
% Author(s): SDG
% Created: 05-Dec-2007
% $Revision: 1 $  $Date: 2007/11/27 14:06:02 $
Theta=zeros(size(Est));
for i=1:size(Est,1)
    for j=1:size(Est,2)
%         if i==310 & j==1
%             disp('Stop')
%         end
        if Nth(i,j)==0 & Est(i,j)==0
            Theta(i,j)=0;
        elseif Nth(i,j)>=0 & Est(i,j)>=0
            if Nth(i,j)==0
                Theta(i,j)=90;
            else
                Theta(i,j)=rad2degFloat(atan(abs(Est(i,j))/abs(Nth(i,j))));
            end
        elseif Nth(i,j)<=0 & Est(i,j)>=0
            if Est(i,j)==0
                Theta(i,j)=180;
            else
                Theta(i,j)=90+rad2degFloat(atan(abs(Nth(i,j))/abs(Est(i,j))));
            end
        elseif Nth(i,j)<=0 & Est(i,j)<=0
            if Nth(i,j)==0
                Theta(i,j)=270;
            else
                Theta(i,j)=180+rad2degFloat(atan(abs(Est(i,j))/abs(Nth(i,j))));
            end
        elseif Nth(i,j)>=0 & Est(i,j)<=0
            if Est(i,j)==0
                Theta(i,j)==0;
            else
                Theta(i,j)=270+rad2degFloat(atan(abs(Nth(i,j))/abs(Est(i,j))));
            end
        end
    end
end
Speed=sqrt(Est.^2+Nth.^2);
return;

function deg=rad2degFloat(rad)

%converts angle in degrees to angle in radians

faktor = pi/180;
deg=rad/faktor;

return;



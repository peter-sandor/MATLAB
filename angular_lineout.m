function histo = angular_lineout(pic,center,Nsectors)
%reads in a monochromatic picture, and its center and calculates its
%velocity distribution sorts pixels into sectors
%Ver3 distinguishes between left and right of image
%in other words, the angle goes from 0 to 2*pi
%angle 0 is defined to correspond to the positive y-axis, angle increases
%counterclockwise
%format for center: [x y] = [column row]
[dimy,dimx,tmp] = size(pic);
tmpx = ((1 : dimx) - center(1)).^2;
tmpy = ((1 : dimy) - center(2)).^2;
Rmax = ceil(max([sqrt(tmpx(1)+tmpy(1)), sqrt(tmpx(1)+tmpy(dimy)), sqrt(tmpx(dimx)+tmpy(1)), sqrt(tmpx(dimx)+tmpy(dimy))]));
histo = zeros([Nsectors Rmax]);
A = histo;
%sector vector
dS = 2*pi/Nsectors;
for indx = 1 : dimx
    for indy = 1 : dimy
        %distance
        R =  floor(sqrt(tmpx(indx) + tmpy(indy)))+1;
        
        %angle and sector for pixel, borders count to both sectors
        alpha = angle(-i*(indy-center(2))+(indx-center(1)))-pi/2;
%         alpha = mod(alpha,2*pi);
        if alpha < 0
            alpha = 2*pi+alpha;
        end
        indS = max(1,ceil(alpha/dS));
        if round(R) == R
            histo(indS,R) = histo(indS,R) + pic(indy,indx);
            A(indS,R) = A(indS,R) + 1;
        else
            R1 = floor(R);
            R2 = ceil(R);
            histo(indS,R1) = histo(indS,R1) + (R2-R)*pic(indy,indx);
            A(indS,R1) = A(indS,R1) + (R2-R);
            histo(indS,R2) = histo(indS,R2) + (R-R1)*pic(indy,indx);
            A(indS,R2) = A(indS,R2) + (R-R1);   
        end;            
    end;
end;
for indS = 1 : Nsectors
    A(indS,:) = (A(indS,:)==0) + A(indS,:);
    histo(indS,:) = histo(indS,:) ./A(indS,:) .*(1:Rmax) ;
end
end
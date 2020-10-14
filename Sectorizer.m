function Sectorizer(A)
%---------------------------------------------------------------------------------------------
%Created by Sean Crowe
%
%The purpose of this program is to break a picture up into radial and
%angular sectors giving the total yield per sector
%
%
%---------------------------------------------------------------------------------------------

% prompt='please enter the name of your picture: ';
% str = input(prompt,'s');
% 
% A=imread(str);
S=size(A);
X_end=S(1);
Y_end=S(2);
Yield=0.0;

for i=1:1:X_end
    
    for j=1:1:Y_end
    
    A_gray(i,j)=(sum(double(A(i,j,:)).^2))^(1.0/3.0);
    end
end%GrayScales RGB image( remove in final Copy of program, data is already grayscale)

A_gray=uint8(A_gray);


for i=1:1:X_end
    
    for j=1:1:Y_end
        
        Yield=Yield+double(A_gray(i,j));
    end
end%Determines total Yield
disp('-------------------')
disp('The total yield is')
disp(Yield),disp('units')
disp('-------------------')

R_Sectors=input('Please enter the number of radial sectors you would like:  ');
A_Sectors=input('Please enter the number of Angular sectors you would like:  ');


big=((A_Sectors+R_Sectors)/20.0)+3.0; % A good size for the picture to reduce error



for i=1:1:X_end
    
    for j=1:1:Y_end
    
        for n=1:1:big
            
            for m=1:1:big
                
                  A_big(big*i+m,big*j+n)=A_gray(i,j);
                  
            end
            
        end
    end
end%Enlarges Picture


image(A_big)
    


d_theta=(pi/2)/A_Sectors;

    if(X_end<Y_end)
    d_r=(X_end*big)/R_Sectors;
    r=(X_end*big);
    else
    d_r=(Y_end*big)/R_Sectors;
    r=Y_end*big;
    end
    
    info=zeros(R_Sectors,A_Sectors);
    
   for rad=0.0:d_r:(r-d_r)
      
      for theta=0.0:d_theta:(pi/2-d_theta)
          
          X_low=floor(cos(theta+d_theta)*rad);
          X_high=ceil(cos(theta)*(rad+d_r));
          Y_high=ceil(sin(theta+d_theta)*(rad+d_r));
          Y_low=floor(sin(theta)*rad);
          
          for i=X_low:1:X_high-1
              
              for j=Y_low:1:Y_high-1
                  
                  if((sqrt(i^2+j^2)>rad)&&((sqrt(i^2+j^2)<(rad+d_r))))
                      
                      if((atan(j/i)>theta)&&((theta+d_theta)>atan(j/i)))
                          
   info(uint8((rad/d_r)+1),uint8(theta/d_theta)+1)=info(uint8((rad/d_r)+1),uint8(theta/d_theta)+1)+(double(A_big(i+1,j+1))/(double(big)^2));
                    
   
                      end
                  end
              end
          end
          
      end
   end%Goes to each sector, determines yield in that sector.
          
          
   
 surf(info)

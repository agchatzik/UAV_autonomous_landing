close all
clc
clear

%{
video = VideoWriter('coke.avi'); %create the video object
open(video); %open the file for writing

N = 291;

for ii=1:N %where N is the number of images
  if(ii<10)  
     filename = [sprintf('000%d',ii) '.jpg'];
  elseif(ii<100)
     filename = [sprintf('00%d',ii) '.jpg'];   
  else
     filename = [sprintf('0%d',ii) '.jpg'];   
  end
  I = imread(filename); %read the next image
  writeVideo(video,I); %write the image to file
end
close(video); %close the file

%}


%video writer
export_video = VideoWriter('video3.avi'); 
open(export_video); %open the file for writing

%initial position of rectangle model
model_center = [320 200];

%define rectangle's size
model_size = [50 80];

%coefficient for the noise covariance
C=2; %bigger C => higher dispersion of particles

%mean and variance for the noise model 
M = [0 0]';
V = C*[2 0.5; 0.5 2];

%number of particles
N = 50; 

%define constant value for likelihood
a = 3;

%initialize particles
particles = floor( mvnrnd(M,V,N) + repmat([model_center(1) model_center(2)], N, 1));

% initialize particle weights to the same value
w = ones(1,N)/N; 

%Read a video
v = VideoReader('Coke.avi');

%Extract Frames from video
video = read(v);
[k,l,n,numFrames]= size(video);

R = (sqrt((model_size(1)^2) + (model_size(2)^2)))/2;
fi = atan((model_size(1)/model_size(2)));  %rad

for frame_count = 1:numFrames-1
    
     fr = video(:,:,:,frame_count);
   
     % downsize to half
     im = im2double(rgb2gray(fr));

     figure;
     out_image = fr;
     imshow(out_image);
         
     edge_image = edge(im,'canny');

     %compute distance trannsform from edges
     dist = bwdist(edge_image);
     
     hold on;
     for i = 1 : N
         count = 0 ;
         dist_sum = 0;
         
         particle = particles(i,:);
         
         %for each particle we create a rectangle with fixed size and
         %center equal to the particles coordinates
        
         rowFrom = floor(particle(2) - model_size(2)/2);
         colFrom = floor(particle(1) - model_size(1)/2);       
         
         rowTo = floor(particle(2) + model_size(2)/2);
         colTo = floor(particle(1) + model_size(1)/2);
         
         if (colTo > 240)
             colTo = 240;
         end
         
         if (colFrom < 1)
             colFrom = 1;
         end
          
         if (rowTo > 180)
             rowTo = 180;
         end
         
         if (rowFrom < 1)
             rowFrom = 1;
         end
         
         %top edge
         row = floor(particle(2) - (model_size(2)/2));
         if (row < 1)
             row = 1;
         elseif (row >180)
             row = 180;
         end
         
         dist_sum = dist_sum + sum(floor(dist(row, colFrom:colTo)));
         count = count + (colTo - colFrom);
        
         %bottom edge
         row = floor(particle(2) + (model_size(2)/2));
         if (row < 1)
             row = 1;
         elseif (row >180)
             row = 180;
         end
         dist_sum = dist_sum + sum(floor(dist(row,  colFrom:colTo)));
         count = count + (colTo - colFrom);
         
         %right edge
         col = floor(particle(1) - (model_size(1)/2));
         if (col < 1)
             col = 1;
         elseif (col >240)
             col = 240;
         end
         dist_sum = dist_sum + sum(floor(dist(rowFrom:rowTo,col)));
         count = count + (rowTo - rowFrom);
         
         %left edge
         col = floor(particle(1) + (model_size(1)/2));
         if (col < 0)
             col = 0;
         elseif (col > 240)
             col = 240;
         end
         dist_sum = dist_sum + sum(floor(dist(rowFrom:rowTo,col)));
         count = count + (rowTo - rowFrom);
         
         total_score = dist_sum/1000;
         w(i) = exp(-a * (total_score));
                     
     end
      
    %weight normalization
    if sum(w)==0
    break;
    end
    w = w/sum(w);
    
     theta = randi([-1 1],1); 
     for k = 1 : N
         particle = particles(k,:);
       
         x0 =  particle(1) + R*sin(fi - theta);
         x1 =  particle(1) + R*sin(pi - fi - theta);
         x2 =  particle(1) + R*sin(pi + fi - theta);
         x3 =  particle(1) + R*sin(- fi - theta);
         
         y0 =  particle(2) + R*cos(fi - theta);
         y1 =  particle(2) + R*cos(pi - fi - theta);
         y2 =  particle(2) + R*cos(pi + fi - theta);
         y3 =  particle(2) + R*cos(- fi - theta);
            
         plot([x3 x2 x1 x0 x3],[y3 y2 y1 y0 y3],'-k');
     end
     
     %figure;
     %imshow(out_image);
     hold off;
     writeVideo(export_video,mat2gray(out_image)); %write the image to file
  
     
    %do resampling based on the weights
    new = zeros(N,2);
    for c=1:N
        s = performSampling(w);
        M1 = [particles(s,1); particles(s,2)] + M; 
        h = mvnrnd(M,V,1)';
        new(c,:) = M1 + h; %generate observations
    end
    
    particles = new;
end
close(export_video); %close the file
%implay('o_sevenup.mpeg');


%Do sampling given the weight function f
function out = performSampling(f)

x = rand;
acc = 0;
i=1;

while 1
    if i > 50     
        break;
    end
    acc = acc + f(i);
    if acc > x
        break;
    end
    i=i+1;
    
end

out=i;

if i>50
    out = 50;
end

end




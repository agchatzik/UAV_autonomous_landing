close all
clc
clear

%video writer
export_video = VideoWriter('export_video.avi'); 
open(export_video); %open the file for writing

%initial position of rectangle model
model_center = [182 97];

%define rectangle's size
model_size = [75 140];

%coefficient for the noise covariance
C=2; %bigger C => higher dispersion of particles

%mean and variance for the noise model 
M = [0 0]';
V = C*[2 0.5; 0.5 2];

%number of particles
N = 50; 

%define constant value for likelihood
a = 2;

%initialize particles
particles = floor( mvnrnd(M,V,N) + repmat([model_center(1) model_center(2)], N, 1));

% initialize particle weights to the same value
w = ones(1,N)/N; 

%Read a video
v = VideoReader('o_sevenup.mpeg');

%Extract Frames from video
video = read(v);
[k,l,n,numFrames]= size(video);


for frame_count = 1:numFrames-1
    
     fr = video(:,:,:,frame_count);
   
     % downsize to half
     im = im2double(rgb2gray(fr));

     %figure;
     out_image = fr;
     %imshow(out_image);
     %title('C=2 ,N=50');   
    
     edge_image = edge(im,'canny');
     %imshow(edge_image);
   
     %compute distance trannsform from edges
     dist = bwdist(edge_image);
     
     hold on;
     for i = 1 : N
         count = 0 ;
         dist_sum = 0;
         
         particle = particles(i,:);
         
         %for each partivle we create a reactangle with fixed size and
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
         
         %left edge
         col = floor(particle(1) - (model_size(1)/2));
         if (col < 1)
             col = 1;
         elseif (col >240)
             col = 240;
         end
         dist_sum = dist_sum + sum(floor(dist(rowFrom:rowTo,col)));
         count = count + (rowTo - rowFrom);
         
         %right edge
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
    

     for k = 1 : N
         particle = particles(k,:);
         
         out_image = insertShape(out_image,'rectangle',[particle(1)-model_size(1)/2 particle(2)-model_size(2)/2 model_size(1) model_size(2)],'LineWidth',1,'color','black');
         %rectangle('Position',[particle(1)-model_size(1)/2 particle(2)-model_size(2)/2 model_size(1) model_size(2)]);
    
     end
     
     %figure;
     %imshow(out_image);
     %title('C=2 ,N=70'); 
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
ii=1;

while 1
    if ii > 50     
        break;
    end
    acc = acc + f(ii);
    if acc > x
        break;
    end
    ii=ii+1;
    
end

out=ii;

if ii>50
    out = 50;
end

end




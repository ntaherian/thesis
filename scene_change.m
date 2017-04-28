function distance = scene_change( frame1,frame2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

size_frame = size(frame1);
frame_num_param = size_frame(1) * size_frame(2) * size_frame(3);

hist_f1 = imhist(frame1);
hist_f1 = hist_f1/frame_num_param;
hist_f2 = imhist(frame2);
hist_f2 = hist_f2/frame_num_param;

cumulative_f1 = cumsum(hist_f1);
cumulative_f2 = cumsum(hist_f2);

distance = sum(abs(cumulative_f1 - cumulative_f2),1);

end


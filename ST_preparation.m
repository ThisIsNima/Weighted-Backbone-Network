function [Dynamic_Adjacency] = ST_preparation(data_set, w_t)

% Inputs:
%   w_t - the window size of temporal segmentation.
%
%   data_set - the window size of temporal segmentation.
% Outputs:
%  Dynamic_Adjacency - The set of adjacency matrices for dynamic connectivity of the network
%
t=length(data_Set{1,1}); %length of time series
tau= t/w_t; % number of temporal segmentations
num_voxel=size(data_set{1,1},1 ); %number of voxels in region of interest

for f=1:length(data_set)
   x_raw=data_set{1, f};
   k=1;
   
   for i=1:w_t:t-w_t   
      x= x_raw(:,i:i+19);
      a(k).k=x;
      k=k+1;
   end

   clear correlations
   clear Adj

   for i=1:tau
      z= a(i);
      for j=1:num_voxel
      for l=1:num_voxel
         temp_corr=corrcoef(z.k(j,:),z.k(l,:));
         correlations(j,l)=temp_corr(1,2);
      end
      end
      CORRz(i).i=correlations;
      disp('i')
      disp(i)
   end

   for i=1:tau
      Adj{i,1}=CORRz(i).i;
   end
   
   Dynamic_Adjacency{f}=Adj;
end

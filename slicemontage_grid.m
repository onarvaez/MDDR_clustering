function slmo = slicemontage_grid(num_clust)


%---
if num_clust == 2
    
slmo.Nrow = 6;
slmo.left_v = [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
          1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 7   7   8   8  ...
          5.5 5.5 0];
      

slmo.bottom_v = [0.75 4.5 4.5 4.5 4.5 3.5 3.5 3.5 3.5 2.5 2.5 2.5 2.5 4.5 4.5 3.5 3.5 2.5 ...
           1   0   1   0    1   0    1   0   1   0    1   0   1  0  2];
%---
elseif num_clust == 3
slmo.Nrow = 7;

slmo.left_v = [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
          1.5 1.5 1.5 2.5 2.5 2.5 3.5 3.5 3.5 4.5 4.5 4.5 7   7   7  8   8   8 ...
          5.5 5.5 5.5 0];
      

slmo.bottom_v = [3.25 5.5 5.5 5.5 5.5 4.5 4.5 4.5 4.5 3.5 3.5 3.5 3.5 5.5 5.5 4.5 4.5 3.5 ...
          2   1   0  2   1   0   2   1   0    2   1   0   2   1   0   2   1   0  2  1  0  1.75];
%---
elseif num_clust == 4
slmo.Nrow = 8;

slmo.left_v = [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
          1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 4.5 4.5 4.5 4.5 7   7   7  7  8   8   8  8 ...
          5.5 5.5 5.5 5.5 0];
      

slmo.bottom_v = [4.25 6.5 6.5 6.5 6.5 5.5 5.5 5.5 5.5 4.5 4.5 4.5 4.5 6.5 6.5 5.5 5.5 4.5 ...
          3   2   1   0  3   2   1   0  3   2   1   0   3   2   1   0   3   2   1   0   3   2   1   0  3  2  1  0  2.75];
%---
elseif num_clust == 5
slmo.Nrow = 9;
    
slmo.left_v = [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
          1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7  8   8   8  8  8 ...
          5.5 5.5 5.5 5.5 5.5 0];
      

slmo.bottom_v = [5.25 7.5 7.5 7.5 7.5 6.5 6.5 6.5 6.5 5.5 5.5 5.5 5.5 7.5 7.5 6.5 6.5 5.5 ...
         4  3   2   1   0  4  3   2   1   0  4  3   2   1   0  4  3   2   1   0  4  3   2   1   0  4  3   2   1   0  4  3  2  1  0  3.75];

%---
elseif num_clust == 6
slmo.Nrow = 10;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
        1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 4.5 4.5 4.5 4.5 4.5 4.5 7 7 7 7  7  7  8  8 8   8  8  8 ...
        5.5 5.5 5.5 5.5 5.5 5.5 0];
    
slmo.bottom_v = [6.25 8.5 8.5 8.5 8.5 7.5 7.5 7.5 7.5 6.5 6.5 6.5 6.5 8.5 8.5 7.5 7.5 6.5 ...
        5 4  3   2   1   0  5 4  3   2   1   0  5 4  3   2   1   0  5 4  3   2   1   0  5 4  3   2   1   0  5 4  3   2   1   0 ...
        5 4  3  2  1  0  4.75];
%---    
elseif num_clust == 7
slmo.Nrow = 11;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
        1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7 7 7 7 7  7  7  8 ...
        8  8 8   8  8  8 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
    
slmo.bottom_v = [7.25 9.5 9.5 9.5 9.5 8.5 8.5 8.5 8.5 7.5 7.5 7.5 7.5 9.5 9.5 8.5 8.5 7.5 ...
        6 5 4  3   2   1   0  6 5 4  3   2   1   0  6 5 4  3   2   1   0  6 5 4  3   2   1   0 6 5 4  3   2   1   0 6 5 4  3   2   1   0 ...
        6 5 4  3  2  1  0  5.75];

%---
elseif num_clust == 8
slmo.Nrow = 12;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
        1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 ...
        7 7 7 7 7 7  7  7 8 8  8  8 8   8  8  8 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
    
slmo.bottom_v = [8.25 10.5 10.5 10.5 10.5 9.5 9.5 9.5 9.5 8.5 8.5 8.5 8.5 10.5 10.5 9.5 9.5 8.5 ...
        7 6 5 4  3   2   1   0 7 6 5 4  3   2   1   0 7 6 5 4  3   2   1   0 7 6 5 4  3   2   1   0 7 6 5 4  3   2   1   0 7 6 5 4  3   2   1   0 ...
        7 6 5 4  3  2  1  0  6.75];
%---
elseif num_clust == 9
slmo.Nrow = 13;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5  1.5 1.5 1.5  2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7  8   8   8  8  8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [9.25 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 9.5 9.5 9.5 9.5 11.5 11.5 10.5 10.5 9.5 ...
          8   7   6  5  4 3 2 1 0   8   7   6  5  4 3 2 1 0   8   7   6  5  4 3 2 1 0   8   7   6  5  4 3 2 1 0 ...
         8   7   6  5  4 3 2 1 0   8   7   6  5  4 3 2 1 0   8   7   6  5  4 3 2 1 0  7.75];
%---

elseif num_clust == 10
slmo.Nrow = 14;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7 7 8   8   8  8  8 8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [10.25 12.5 12.5 12.5 12.5 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 12.5 12.5 11.5 11.5 10.5 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0  8.75];
     
%---
elseif num_clust == 11
slmo.Nrow = 14;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7 7 8   8   8  8  8 8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [10.25 12.5 12.5 12.5 12.5 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 12.5 12.5 11.5 11.5 10.5 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0  8.75];
     
%---
elseif num_clust == 12
slmo.Nrow = 14;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7 7 8   8   8  8  8 8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [10.25 12.5 12.5 12.5 12.5 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 12.5 12.5 11.5 11.5 10.5 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0  8.75];
  
     
 %---
elseif num_clust == 13
slmo.Nrow = 14;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7 7 8   8   8  8  8 8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [10.25 12.5 12.5 12.5 12.5 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 12.5 12.5 11.5 11.5 10.5 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0  8.75];
     

%---
elseif num_clust == 14
slmo.Nrow = 14;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7 7 8   8   8  8  8 8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [10.25 12.5 12.5 12.5 12.5 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 12.5 12.5 11.5 11.5 10.5 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0  8.75];
     
%---
elseif num_clust == 15
slmo.Nrow = 14;

slmo.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 7   8   7   8   7   ...
         1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 ...
         4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 7   7   7  7  7 7 7 7 7 7 8   8   8  8  8 8 8 8 8 8 ...
         5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 0];
slmo.bottom_v = [10.25 12.5 12.5 12.5 12.5 11.5 11.5 11.5 11.5 10.5 10.5 10.5 10.5 12.5 12.5 11.5 11.5 10.5 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0 ...
         9 8   7   6  5  4 3 2 1 0  9 8   7   6  5  4 3 2 1 0 9  8   7   6  5  4 3 2 1 0  8.75];
     

     

end  

function tetra = calculate_tetrahedron(i,j,k,k_d,j_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% This function returns the 6 tetrahedrons in a cube with first index (i,j,k)
%% It gives the connections between the points in a cube.
%% 
%%  [1,3,8,4;2,1,8,4;2,5,1,4;2,6,5,4;2,8,6,4;2,7,5,6;]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tetra = [k + k_d*(j-1 + j_d*(i-1))    k + k_d*(j + j_d*(i-1))     k+1 + k_d*(j + j_d*(i))      k+1 + k_d*(j + j_d*(i-1)); ...
         k+1 + k_d*(j-1 + j_d*(i-1))  k + k_d*(j-1 + j_d*(i-1))   k+1 + k_d*(j-1 + j_d*(i-1))  k + k_d*(j-1 + j_d*(i-1)); ...
         k+1 + k_d*(j-1 + j_d*(i-1))  k + k_d*(j-1 + j_d*(i))     k + k_d*(j-1 + j_d*(i-1))    k+1 + k_d*(j + j_d*(i-1)); ...
         k+1 + k_d*(j-1 + j_d*(i-1))  k+1 + k_d*(j-1 + j_d*(i))   k + k_d*(j-1 + j_d*(i))      k+1 + k_d*(j + j_d*(i-1)); ...
         k+1 + k_d*(j-1 + j_d*(i-1))  k+1 + k_d*(j + j_d*(i))     k+1 + k_d*(j-1 + j_d*(i))    k+1 + k_d*(j + j_d*(i-1)); ...
         k+1 + k_d*(j-1 + j_d*(i-1))  k + k_d*(j + j_d*(i))       k + k_d*(j-1 + j_d*(i))      k+1 + k_d*(j-1 + j_d*(i))];
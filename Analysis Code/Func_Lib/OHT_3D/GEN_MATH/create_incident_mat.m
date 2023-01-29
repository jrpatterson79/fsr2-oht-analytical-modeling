function [E] = create_incident_mat(num_x, num_y, varargin)
           
%Code for creating an oriented incidence matrix 
%NOTES: 
%-Currently only 2D
%-Currently does not calculate length between two adjacent blocks (assumes
%uniform grid in both directions)
%-Needs to be tested to make sure indexing is right with respect to
%unknowns created by meshgrid and then reshaped
%***Need to think through other methods of doing this (inputs, etc).
%Most user-friendly input would be to take grid and ordered midpoints from
%user (??)

dim = nargin;

if dim == 2
    num_edges = (num_x - 1)*num_y + (num_y - 1)*num_x;
    num_blocks = num_x*num_y;
    E = sparse(num_edges,num_blocks);
    edge_num = 1;
    for j = 1:1:num_x
        for i = 1:1:num_y
            if i < num_y
                first = sub2ind([num_y num_x], i,j);
                second = sub2ind([num_y num_x], i+1,j);
                E(edge_num,first) = 1;
                E(edge_num,second) = -1;
                edge_num = edge_num + 1;
            end
            if j < num_x
                first = sub2ind([num_y num_x], i,j);
                second = sub2ind([num_y num_x], i,j+1);
                E(edge_num,first) = 1;
                E(edge_num,second) = -1;
                edge_num = edge_num + 1;
            end
        end
    end
end
                
    
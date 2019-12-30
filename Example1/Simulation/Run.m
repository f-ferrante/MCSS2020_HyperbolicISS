 
%Run the solver and generates the solution vector [x1, x2] and the grid of
%points (t,z) over which the solution has been computed

[z1,t1, x1, x2]=hyp_static_bc(); 

z=z1(1,:); %extract a vector of samples for the spatial variable z
t=t1(1,:); %extract a vector of samples for the time variable t

Norm=SpatialNorm(x1, x2,z); %Compute the evolution over time of the (spatial) L2 norm of the solution (V1to(t,.),V2to(t,.))

plot(t,Norm,'-b','linewidth', 2); %Plot the (spatial) L2 norm of the solution (V1to(t,.),V2to(t,.))
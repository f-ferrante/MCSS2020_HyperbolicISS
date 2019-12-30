 
%Run the solver and generate the solution vector [V1to, V2to,V3to] and the grid of
%points (t,z) over which the solution has been computed

[x1,t1, V1to, V2to, V3to]=hyp_static_bc(); 

z=x1(1,:);

Norm=SpatialNorm(V1to, V2to,V3to, z); %Compute the evolution over time of the (spatial) L2 norm of the solution (V1to(t,.),V2to(t,.),V3to(t,.))
plot(t1(:,1),Norm,'-k','linewidth', 2); %Plot the evolution over time of the (spatial) L2 norm of the solution (V1to(t,.),V2to(t,.),V3to(t,.))
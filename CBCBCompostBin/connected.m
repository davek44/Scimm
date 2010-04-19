function conn = connected(A)
  % v will serve as a the connectivity from vertex 1
  v = A(1,:);

  % track the size of v in order to assess progress
  vsize = nnz(v);

  notstable = true;
  while(notstable)

    % for each reachable vertex
    [x,y] = find(v);
    for i=1:size(x,2)

      % find reachable vertexes
      [xr,yr] = find(A(y(i),:));
      for j=1:size(xr,2)
        v(yr(j)) = 1;
      end
    end
    
    % changes?
    if(nnz(v) > vsize)
      notstable = true;
    else
      notstable = false;
    end
    vsize = nnz(v);
  end

  conn = (vsize == size(A,1));

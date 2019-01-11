while(1)

    inner_flag = (rand() > 0.5) + 0.0;
    n = 5 + randi(10);
    time_points = randi(10);

    while(1)
        T = zeros(n);
        for i = 2:n
            T(randi(n),i) = 1;
        end
        adj = T + T';
        if (graphisspantree(sparse(adj)) == 1)
            break;
        end
    end

    U = inv(eye(n) - T);
    D =  diag(0.0001 + rand(n,1)); 

    F = rand(n,time_points);
    F_tilde = D*F;
    F_tilde_tilde = D*D*F;
    U_tilde = D*U;
    tiln = U'*F;
    AdjT = T;

    cvx_begin quiet
        cvx_precision best
        variable M(n,time_points)
        minimize (  sum(sum(( D*F - D*U*M ).^2))      )
        subject to
            if (inner_flag == 1) 
                sum(M,1) <= 1
            else
                sum(M,1) == 1
            end
            M >= 0
    cvx_end

    fid = fopen('inputdata.txt','w');

    fprintf(fid, '%d %d\n',n,time_points);
    F_T = F';
    fprintf(fid, '%1.6f ',F_T(1:end-1)); fprintf(fid, '%1.6f\n',F_T(end));
    fprintf(fid, '%1.6f ',D(1:n+1:end-1));fprintf(fid, '%1.6f\n',D(end));
    fprintf(fid, '0\n');
    deg = sum(adj);
    fprintf(fid, '%d ',deg(1:end-1)); fprintf(fid, '%d\n',deg(end));
    for i = 1:n
        neigh = find(adj(:,i))'-1;
        if (length(neigh) > 1)
            fprintf(fid, '%d ',neigh(1:end-1));
        end
        fprintf(fid, '%d\n',neigh(end));
    end
    fprintf(fid, '1\n');

    fclose(fid);

    system(['./a.out inputdata.txt outputdata.txt ', num2str(inner_flag)]);

    fid = fopen('outputdata.txt','r');

    cost_exact = inf;
    M_exact = zeros(n*time_points,1);
    cost_exact = fscanf(fid, '%f\n',1);
    M_exact = fscanf(fid, '%f ',n*time_points);

    M_exact = reshape(M_exact, n,time_points);

    error_sol = max(max(M - M_exact));

    error_obj = abs(cost_exact - cvx_optval);

    disp([n, time_points, inner_flag, error_sol, error_obj]);

    if (error_obj > 0.001 || isnan(error_obj) == 1)
        disp(['NUMERICAL ERROR']);
        break;
    end

end
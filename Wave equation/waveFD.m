% Finite-difference for wave equation with specified initial conditions
% (on retangular domain)
%
% @param: space interval [x_start, x_end],
%         time interval [t_start, t_end],
%         number of space steps M,
%         number of time steps n


function w = waveFD(x_start,x_end,t_start,t_end,M,n)
  % Initial info of the wave PDE
  f = @(x) log(1 + x.^1);
  f_t = @(x) 1./(1 + x.^1);
  l = @(t) log(1 + t.^1);
  r = @(t) log(2 + t.^1);
  D = 1;

  % Setting up step sizes and alpha = D*k/h
  h = (x_end - x_start)/M;
  k = (t_end - t_start)/n;
  m = M-1;
  alpha = (D*k)/h;

  % Input parameters check
  if M < 4 || n < 4 || x_end - x_start < 0 || t_end - t_start < 0
    printf("Invalid parameters inputted\n");
    return;
  end

  % CFL condition check
  if alpha > 1
    printf("Fail to meet the CFL condition for stability with the inputted steps\n");
    return;
  end

  % Constructing finite-difference matrix
  A = diag(2-2*(alpha^2)*ones(m,1)) ...
    + diag((alpha^2)*ones(m-1,1),1) ...
    + diag((alpha^2)*ones(m-1,1),-1);

  % Initialize the solution
  w(:,1) = (1/2)*A*transpose(f(x_start + (1:m)*h)) + k*transpose(f_t(x_start + (1:m)*h));
  w(:,2) = A*w(:,1) - transpose(f(x_start + (1:m)*h));

  % Iterations
  for j = 2:n
    w(:,j+1) = A*w(:,j) - w(:,j-1);
  end

  % Boundary attaching
  w = [alpha^2*l((0:n)*k + t_start) ; w ; alpha^2*r((0:n)*k + t_start)];
  w(1,1) = w(1,1)/2;
  w(m+2, 1) = w(m+2,1)/2;

  % Plotting
  x = (0:m+1)*h + x_start;
  t = (0:n)*k + t_start;
  mesh(x,t,transpose(w))
  xlabel('x');
  ylabel("t");
  zlabel("u(x,t)")
  view(60,30);
  axis([x_start, x_end, t_start, t_end, min(w(:))-0.25, max(w(:))+0.25])

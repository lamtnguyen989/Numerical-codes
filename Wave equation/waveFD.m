% Finite-difference for wave equation with specified initial conditions
% (on retangular domain)
%
% @param: space interval [x_0, x_end],
%         time interval [t_0, t_end],
%         number of space steps m,
%         number of time steps n

function w = waveFD(x_0,x_end,t_0,t_end,m,n)
  % Initial info of the wave PDE
  f = @(x) 0*x;
  l = @(t) 0*t;
  r = @(t) 0*t;
  f_t = @(x) 2*pi*sin(pi*x);
  D = 2;

  % Setting up step sizes and alpha = D*k/h
  h = (x_end-x_0)/m;
  k = (t_end-t_0)/n;
  alpha = (D*k)/h;

  % CFL condition check
  if alpha > 1
    printf("Fail to meet the CFL condition for stability\n")
    return;
  end

  % Constructing finite-difference matrix
  A = diag(2-2*(alpha^2)*ones(m,1)) + diag((alpha^2)*ones(m-1,1),1) + diag((alpha^2)*ones(m-1,1),-1);

  % Initialize the solution
  w(:,1) = transpose(f(x_0+(1:m)*h));
  w(:,2) = (1/2)*A*w(:,1) + k*transpose(f_t(x_0 + (1:m)*h)) + (1/2)*(alpha^2)*[l(t_0); zeros(m-2,1); r(t_0)];

  % Iterations
  for j = 2:n
    w(:,j+1) = A*w(:,j) - w(:,j-1) + (alpha^2)*[l(t_0 + k*j); zeros(m-2,1); r(t_0 + k*j)];
  end

  % Plotting
  x = (0:m-1)*h + ones(1,m)*x_0;
  t = (0:n)*k + ones(1,n+1)*t_0;
  mesh(x,t,w')
  xlabel("x");
  ylabel("t");
  zlabel("u(x,t)")
  view(60,30);
  axis([x_0, x_end, t_0, t_end, min(w(:))-0.25, max(w(:))+0.25])

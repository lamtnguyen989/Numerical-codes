% Simpsion's rule calculation
function res = simpson(f,start,endpt,m)
  % Step size
  h = (endpt - start)/(2*m);

  % Initializing the result of the computation
  res = f(start) + f(endpt);

  % Making initial jump
  x = start + h;

  % Sum calculation
  for i = 1 : 2*m-1
    if mod(i,2) == 1
      res = res + 4*f(x);
    else
      res = res + 2*f(x);
    end
    % Iteratively adding step-size to gain new inputs to calculate the sum
    x = x + h;
  end
  % Tie-up loose-end
  res = (h/3)*res;

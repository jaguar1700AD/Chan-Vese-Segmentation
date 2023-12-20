pkg load image;

file = "apple.jpeg";

% Set constants
neta = 10^(-8);
mu = 0.2;
nu = 0;
lambda1 = 1; lamnda2 = 1;
dt = 0.5;
tol = 0.6*10^(-2);
global epsilon = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function display_image(f)
  figure(1); clf;
  if (length(size(f)) == 2)
    colormap(gray);
  endif
  imagesc(f);
  %imshow(f);
  caxis([0 1])
endfunction

function result = dirac_delta(A)
  global epsilon;
  result = epsilon ./ (pi * (epsilon^2 + A.^2));
endfunction

function result = heaviside(A)
  global epsilon;
  result = (1 + (2 * atan(A / epsilon) / pi)) / 2;
endfunction

function result = integrate(A)
  result = sum(sum(A));
endfunction

function result = average(f, density)
  result = integrate(f .* density) / integrate(density);
endfunction

%Extend f to outside the image boundary using 0 Neumann Boundary Condition
function result = extend(f)
  A = [f(:,1, :), f, f(:,end, :)];
  result = [A(1,:, :); A; A(end, :,:)];
endfunction

function result = restrict(f)
  result = f(2:end-1, 2:end-1,:);
endfunction

%%% Functions for various differences. Coded in such a way that the output
%% is 0 whenever a difference doesn't make sense

function result = row_forward_diff(f)
  result = [f(2:end,:,:); f(end,:,:)] - f;
endfunction

function result = row_backward_diff(f)
  result = f - [f(1,:,:); f(1:end-1,:,:)];
endfunction

function result = row_central_diff(f)
  result = ([f(2:end, :,:); f(end-1,:,:)] - [f(2,:,:); f(1:end-1,:,:)]) / 2;
endfunction

function result = col_forward_diff(f)
  result = [f(:, 2:end,:), f(:, end,:)] - f;
endfunction

function result = col_backward_diff(f)
  result = f - [f(:, 1,:), f(:, 1:end-1,:)];
endfunction

function result = col_central_diff(f)
  result = ([f(:, 2:end,:), f(:, end-1,:)] - [f(:, 2,:), f(:, 1:end-1,:)]) / 2;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read an image from a file
f = imread(file);
info = imfinfo(file);

% convert image to double and scale to [0,1]
f = double(f);
f = f / 2^(info.BitDepth);

[n1,n2, n3] = size(f)
LargestPixel = norm(reshape(f, n1*n2*n3, 1), inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize phi
x = 1:n2; y = (1:n1)';
[xx, yy] = meshgrid(x, y);
phi = sin(pi * 10* x / n2) .* sin(pi * 10* y / n1);

display_image(f);

step = 1;
while (1)

  phi_extended = extend(phi); % Extend using 0 neumann bdd condns

  phi_i_plus = row_forward_diff(phi_extended);
  phi_i_minus = row_backward_diff(phi_extended);
  phi_i_zero = row_central_diff(phi_extended);

  phi_j_plus = col_forward_diff(phi_extended);
  phi_j_minus = col_backward_diff(phi_extended);
  phi_j_zero = col_central_diff(phi_extended);

  A = mu ./ ((neta^2 + phi_i_plus.^2 + phi_j_zero.^2).^(0.5));
  B = mu ./ ((neta^2 + phi_i_zero.^2 + phi_j_plus.^2).^(0.5));

  d = heaviside(phi);
  c1 = average(f, d);
  c2 = average(f, 1 - d);

  % Do all difference computations using phi_extended to enforce 0 Neumann bdd condns
  Dphi = restrict(row_backward_diff(A .* phi_i_plus) + col_backward_diff(B .* phi_j_plus));
  Dphi = Dphi - nu - lambda1 * sum((f - c1).^2,3) + lamnda2 * sum((f - c2).^2, 3);
  Dphi = Dphi .* dirac_delta(phi);

  phi_new = phi + Dphi * dt;

  if (rem(step, 100) == 0)
    step
    display_image(f);
    hold on;
    contour(x, y, phi_new, [0,0], 'r');
    drawnow
    pause(0.001);
    if (norm(reshape(phi_new - phi, n1*n2, 1), 2) / sqrt(n1*n2) < tol)
      printf("Converged\n");
      break;
    endif
  endif

  phi = phi_new;
  step += 1;

endwhile


%% Output original and result side-by-side
result = [f, f];
imwrite(result, 'result.png')

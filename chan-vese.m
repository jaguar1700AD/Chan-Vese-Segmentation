%% Read an image from a file
f = imread('testpat_blur4.png');

% Set constants
neta = 10^(-8);
mu = 0.1;
global epsilon = 1;

% convert image to double and scale to [0,1]
f = double(f) / 255;
[n1,n2] = size(f);

function display_image(f)
  figure(1); clf;
  %pcolor(u); % try this one too
  imagesc(f);
  caxis([0 1])
  colormap(gray)
  axis equal, axis tight
endfunction

function result = dirac_delta(A)
  global epsilon;
  result = epsilon ./ (pi * (epsilon^2 + A.^2));
endfunction

function result = heaviside(A)
  global epsilon;
  result = (1 + (2 * atan(A / epsilon) / pi)) / 2;
endfunction

%Extend f to outside the image boundary using 0 Neumann Boundary Condition
function result = extend(f)
  A = [f(:,1), f, f(:,end)];
  result = [A(1,:); A; A(end, :)];
endfunction

function result = restrict(f)
  result = f(2:end-1, 2:end-1);
endfunction

%%% Functions for various differences. Coded in such a way that the output
%% is 0 whenever a difference doesn't make sense

function result = row_forward_diff(f)
  result = [f(2:end,:); f(end,:)] - f;
endfunction

function result = row_backward_diff(f)
  result = f - [f(1,:); f(1:end-1,:)];
endfunction

function result = row_central_diff(f)
  result = ([f(2:end, :); f(end-1,:)] - [f(2,:); f(1:end-1,:)]) / 2;
endfunction

function result = col_forward_diff(f)
  result = [f(:, 2:end), f(:, end)] - f;
endfunction

function result = col_backward_diff(f)
  result = f - [f(:, 1), f(:, 1:end-1)];
endfunction

function result = col_central_diff(f)
  result = ([f(:, 2:end), f(:, end-1)] - [f(:, 2), f(:, 1:end-1)]) / 2;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize phi
[x, y] = meshgrid(1:n2, 1:n1);
phi = sin(pi * x / 5) * sin(pi * y / 5);
phi = extend(phi);

for step = 1:100

  phi_i_plus = row_forward_diff(phi);
  phi_i_minus = row_backward_diff(phi);
  phi_i_zero = row_central_diff(phi);

  phi_j_plus = col_forward_diff(phi);
  phi_j_minus = col_backward_diff(phi);
  phi_j_zero = col_central_diff(phi);

  A = mu ./ (neta^2 + phi_i_plus.^2 + phi_j_zero.^2)^(0.5);
  B = mu ./ (neta^2 + phi_i_zero.^2 + phi_j_plus.^2)^(0.5);

endfor


%% Output original and result side-by-side
result = [f, f];
imwrite(result, 'result.png')

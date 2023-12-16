%% Read an image from a file
f = imread('testpat_blur4.png');

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
  epsilon = 1;
  result = epsilon ./ (pi * (epsilon^2 + A.^2));
endfunction

%Extend f to outside the image boundary using 0 Neumann Boundary Condition
function result = extend(f)
  A = [f(:,1), f, f(:,end)];
  result = [A(1,:); A; A(end, :)];
endfunction

function result = row_forward_diff(f)
  result = [f(2:end,:); f(end,:)] - f;
endfunction

function result = row_backward_diff(f)
  result = f - [f(1,:); f(1:end-1,:)];
endfunction

function result = col_forward_diff(f)
  result = [f(:, 2:end), f(:, end)] - f;
endfunction

function result = col_backward_diff(f)
  result = f - [f(:, 1), f(:, 1:end-1)];
endfunction

%% Output original and result side-by-side
result = [f, f];
imwrite(result, 'result.png')

function make_noisy(fname, extension, thresh)
  file = [fname "." extension];
  f = imread(file);
  info = imfinfo(file);
  f = double(f);
  f = f / 2^(info.BitDepth);
  [n1,n2, n3] = size(f);

  % add noise
  u = f + thresh*(rand(size(f)) - 1/2);
  display_image(u);

  imwrite(u, [fname "_" num2str(thresh) "_noisy." extension]);
endfunction

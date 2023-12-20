% Uses canny edge detection from ocatave's image package

function canny(file, thresh)
  f = imread(file);
  canny_edges = edge(f, "Canny", thresh);
  display_image(canny_edges);
endfunction



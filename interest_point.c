#include <stdlib.h>
#include <stdio.h>
#include "Util.h"
#include <math.h>
#include <string.h>

int maxval = 255;
float *sobel_x, *sobel_y, *binomial_3;

float gradient_1d(gray *map, int i, int j, int rows, int cols, float *mask) {
  if (i - 1 < 0 || i + 1 >= rows || j - 1 < 0 || j + 1 >= cols) {
    return 0;
  }

  float val = 0.0f;
  for (int i1 = i - 1; i1 <= i + 1; i1++) {
    for (int i2 = j - 1; i2 <= j + 1; i2++) {
      int mask_index = (i1 - i + 1) * (2 * 1 + 1) + (i2 - j + 1);
      int img_index = i1 * cols + i2;
      val += (float) map[img_index] * mask[mask_index];
    }
  }
  return val;
}

float gradient_x2(gray *map, int i, int j, int rows, int cols) {
  float gradient_x = gradient_1d(map, i, j, rows, cols, sobel_x);
  return gradient_x * gradient_x;
}

float gradient_y2(gray *map, int i, int j, int rows, int cols) {
  float gradient_y = gradient_1d(map, i, j, rows, cols, sobel_y);
  return gradient_y * gradient_y;
}

float gradient_xy(gray *map, int i, int j, int rows, int cols) {
  float gradient_x = gradient_1d(map, i, j, rows, cols, sobel_x);
  float gradient_y = gradient_1d(map, i, j, rows, cols, sobel_y);
  return gradient_x * gradient_y;
}

gray *float_to_image(float *gradient, int rows, int cols) {
  float min_mag = 0;
  float max_mag = 0;
  for (int i=0; i<rows * cols; i++) {
    if (gradient[i] > max_mag) {
      max_mag = gradient[i];
    }
    if (gradient[i] < min_mag) {
      min_mag = gradient[i];
    }
  }

  gray *out = malloc(rows * cols * sizeof(gray));
  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      out[i * cols + j] = (gray) ((gradient[i * cols + j] - min_mag) / (max_mag - min_mag) * maxval);
    }
  }
  return out;
}

float filter_gradient_pixel(float *gradient, int i, int j, int rows, int cols, float *mask) {
  int index = i * cols + j;
  if (i - 1 < 0 || i + 1 >= rows || j - 1 < 0 || j + 1 >= cols) {
    return gradient[index];
  }

  float val = 0.0f;
  for (int i1 = i - 1; i1 <= i + 1; i1++) {
    for (int i2 = j - 1; i2 <= j + 1; i2++) {
      int mask_index = (i1 - i + 1) * (2 * 1 + 1) + (i2 - j + 1);
      int img_index = i1 * cols + i2;
      val += (float) gradient[img_index] * mask[mask_index];
    }
  }
  return val;
}

float *filter_gradient(float *gradient, int rows, int cols, float *mask) {
  float *out = malloc(rows * cols * sizeof(float));

  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      out[i * cols + j] = filter_gradient_pixel(gradient, i, j, rows, cols, mask);
    }
  }

  return out;
}

float *harris(float *grad_x2, float *grad_y2, float *grad_xy, float alpha, int rows, int cols) {
  float *out = malloc(rows * cols * sizeof(float));
  int index;
  float det;
  float trace;

  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      index = i * cols + j;
      det = grad_x2[index] * grad_y2[index] - grad_xy[index] * grad_xy[index];
      trace = grad_x2[index] + grad_y2[index];

      out[index] = det - alpha * trace * trace;
    }
  }

  return out;
}


float *sobel_mask_x() {
  printf("Sobel x\n");
  float *mask = malloc(sizeof(float) * 9);
  mask[0] = -1.0f / 4;
  mask[1] = 0;
  mask[2] = 1.0f / 4;
  mask[3] = -2.0f / 4;
  mask[4] = 0;
  mask[5] = 2.0f / 4;
  mask[6] = -1.0f / 4;
  mask[7] = 0;
  mask[8] = 1.0f / 4;
  return mask;
}

float *sobel_mask_y() {
  printf("Sobel y\n");
  float *mask = malloc(sizeof(float) * 9);
  mask[0] = -1.0f / 4;
  mask[1] = -2.0f / 4;
  mask[2] = -1.0f / 4;
  mask[3] = 0;
  mask[4] = 0;
  mask[5] = 0;
  mask[6] = 1.0f / 4;
  mask[7] = 2.0f / 4;
  mask[8] = 1.0f / 4;
  return mask;
}

float *binomial_3_mask() {
  printf("Binomial 3x3\n");
  float *mask = malloc(sizeof(float) * 9);
  mask[0] = 1.0f / 16;
  mask[1] = 2.0f / 16;
  mask[2] = 1.0f / 16;
  mask[3] = 2.0f / 16;
  mask[4] = 4.0f / 16;
  mask[5] = 2.0f / 16;
  mask[6] = 1.0f / 16;
  mask[7] = 2.0f / 16;
  mask[8] = 1.0f / 16;
  return mask;
}

float *maxima(float *array, float *sorted_array, int nb_maxima, int rows, int cols) {
  int size = rows * cols;
  float *out = malloc(size * sizeof(float));
  for (int i = 0; i < size; i++) {
    if (array[i] < sorted_array[size - nb_maxima]) {
      out[i] = 0;
    } else {
      out[i] = array[i];
    }
  }
  return out;
}

float *suppress_non_local_max(float *array, int rows, int cols) {
  float *out = malloc(rows * cols * sizeof(float));
  int index;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      index = i * cols + j;
      if (array[index] > array[(i - 1) * cols + j - 1]
        && array[index] > array[(i - 1) * cols + j]
        && array[index] > array[(i - 1) * cols + j + 1]
        && array[index] > array[i * cols + j - 1]
        && array[index] > array[i * cols + j + 1]
        && array[index] > array[(i + 1) * cols + j - 1]
        && array[index] > array[(i + 1) * cols + j]
        && array[index] > array[(i + 1) * cols + j + 1]
      ) {
        out[index] = array[index];
      } else {
        out[index] = 0;
      }
    }
  }
  return out;
}

int main(int argc, char* argv[]) {
  FILE* ifp, *ofp;
  gray* graymap;
  int ich1, ich2, rows, cols, pgmraw, nb_maxima;
  int i, j;
  float harris_alpha;

  sobel_x = sobel_mask_x();
  sobel_y = sobel_mask_y();
  binomial_3 = binomial_3_mask();


  /* Arguments */
  if ( argc != 5 ){
    printf("\nUsage: %s [file-in] [file-out] [harris-alpha] [nb_maxima]\n\n", argv[0]);
    exit(0);
  }

  /* Opening */
  ifp = fopen(argv[1],"r");
  if (ifp == NULL) {
    printf("error in opening file %s\n", argv[1]);
    exit(1);
  }

  ofp = fopen(argv[2], "w");
  if (ifp == NULL) {
    printf("error in opening file %s\n", argv[2]);
    exit(1);
  }

  harris_alpha = atof(argv[3]);
  nb_maxima = atoi(argv[4]);

  /*  Magic number reading */
  ich1 = getc( ifp );
  if ( ich1 == EOF )
  pm_erreur( "EOF / read error / magic number" );
  ich2 = getc( ifp );
  if ( ich2 == EOF )
  pm_erreur( "EOF /read error / magic number" );
  if(ich2 != '2' && ich2 != '5')
  pm_erreur(" wrong file type ");
  else
  if(ich2 == '2')
  pgmraw = 0;
  else pgmraw = 1;

  /* Reading image dimensions */
  cols = pm_getint( ifp );
  rows = pm_getint( ifp );
  maxval = pm_getint( ifp );

  /* Memory allocation  */
  graymap = (gray *) malloc(cols * rows * sizeof(gray));

  /* Reading */
  for(i=0; i < rows; i++)
  for(j=0; j < cols ; j++)
  if(pgmraw)
  graymap[i * cols + j] = pm_getrawbyte(ifp) ;
  else
  graymap[i * cols + j] = pm_getint(ifp);

  /* Writing */
  fprintf(ofp, "P2\n");

  fprintf(ofp, "%d %d \n", cols, rows);
  fprintf(ofp, "%d\n",maxval);

  float *grad_x2_image = malloc(cols * rows * sizeof(float));
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      grad_x2_image[i * cols + j] = gradient_x2(graymap, i, j, rows, cols);
    }
  }

  float *grad_y2_image = malloc(cols * rows * sizeof(float));
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      grad_y2_image[i * cols + j] = gradient_y2(graymap, i, j, rows, cols);
    }
  }

  float *grad_xy_image = malloc(cols * rows * sizeof(float));
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      grad_xy_image[i * cols + j] = gradient_xy(graymap, i, j, rows, cols);
    }
  }

  float *grad_x2_image_filtered = filter_gradient(grad_x2_image, rows, cols, binomial_3);
  float *grad_y2_image_filtered = filter_gradient(grad_y2_image, rows, cols, binomial_3);
  float *grad_xy_image_filtered = filter_gradient(grad_xy_image, rows, cols, binomial_3);

  // float *harris_image = harris(grad_x2_image, grad_y2_image, grad_xy_image, harris_alpha, rows, cols);
  float *harris_image = harris(grad_x2_image_filtered, grad_y2_image_filtered, grad_xy_image_filtered, harris_alpha, rows, cols);

  float *sorted_harris = malloc(cols * rows * sizeof(float));
  memcpy(sorted_harris, harris_image, cols * rows * sizeof(float));
  qsort(sorted_harris, cols * rows, sizeof(float), cmp_float);

  float *non_local_harris = suppress_non_local_max(harris_image, rows, cols);
  float *sorted_non_local_harris = malloc(cols * rows * sizeof(float));
  memcpy(sorted_non_local_harris, non_local_harris, cols * rows * sizeof(float));
  qsort(sorted_non_local_harris, cols * rows, sizeof(float), cmp_float);

  float *harris_maxima = maxima(harris_image, sorted_harris, nb_maxima, rows, cols);
  float *harris_lmaxima = maxima(non_local_harris, sorted_non_local_harris, nb_maxima, rows, cols);
  // gray *out = float_to_image(grad_x2_image, rows, cols);
  // gray *out = float_to_image(grad_y2_image, rows, cols);
  // gray *out = float_to_image(grad_xy_image, rows, cols);
  // gray *out = float_to_image(grad_x2_image_filtered, rows, cols);
  // gray *out = float_to_image(grad_y2_image_filtered, rows, cols);
  // gray *out = float_to_image(grad_xy_image_filtered, rows, cols);
  // gray *out = float_to_image(harris_image, rows, cols);
  // gray *out = float_to_image(harris_maxima, rows, cols);
  gray *out = float_to_image(harris_lmaxima, rows, cols);

  for (i=0; i<rows; i++) {
    for(j=0; j<cols; j++) {
      fprintf(ofp, "%d ", out[i * cols + j]);
    }
    fprintf(ofp, "\n");
  }


  /* Closing */
  fclose(ifp);
  return 0;
}

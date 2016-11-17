#include <stdlib.h>
#include <stdio.h>
#include "Util.h"
#include <math.h>

int maxval = 255;
float *sobel_x, *sobel_y;

float gradient_1d(gray *map, int i, int j, int rows, int cols, float *mask) {
  if (i - 1 < 0 || i + 1 > rows || j - 1 < 0 || j + 1 > cols) {
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

gray *filter(float *gradient, int rows, int cols) {
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


int main(int argc, char* argv[]) {
  FILE* ifp, *ofp;
  gray* graymap;
  int ich1, ich2, rows, cols, pgmraw;
  int i, j;

  sobel_x = sobel_mask_x();
  sobel_y = sobel_mask_y();


  /* Arguments */
  if ( argc != 3 ){
    printf("\nUsage: %s file \n\n", argv[0]);
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

  gray *out = filter(grad_xy_image, rows, cols);

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

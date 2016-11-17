#include <stdlib.h>
#include <stdio.h>
#include "Util.h"
#include <math.h>

float *sobel_x, *sobel_y;

float filter_pixel(gray *map, int i, int j, int rows, int cols, float *mask, int neighbor_size) {
  if (i - neighbor_size < 0 || i + neighbor_size > rows || j - neighbor_size < 0 || j + neighbor_size > cols) {
    return 0;
  }

  float val = 0.0f;
  for (int i1 = i - neighbor_size; i1 <= i + neighbor_size; i1++) {
    for (int i2 = j - neighbor_size; i2 <= j + neighbor_size; i2++) {
      int mask_index = (i1 - i + neighbor_size) * (2 * neighbor_size + 1) + (i2 - j + neighbor_size);
      int img_index = i1 * cols + i2;
      val += (float) map[img_index] * mask[mask_index];
    }
  }
  return val;
}

gray *filter(gray *map, int rows, int cols, float *mask, int neighbor_size) {
  gray base = 127;
  float *gradient = malloc(rows * cols * sizeof(float));

  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      gradient[i * cols + j] = filter_pixel(map, i, j, rows, cols, mask, neighbor_size);
    }
  }

  float max_mag = 0;
  for (int i=0; i<rows * cols; i++) {
    if (gradient[i] > 0 && gradient[i] > max_mag) {
      max_mag = gradient[i];
    } else if (gradient[i] < 0 && -gradient[i] > max_mag) {
      max_mag = -gradient[i];
    }
  }

  gray *out = malloc(rows * cols * sizeof(gray));
  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      out[i * cols + j] = (gray) ((float) base + gradient[i * cols + j] / max_mag * base);
    }
  }
  return out;
}

float filter_pixel_mag(gray *map, int i, int j, int rows, int cols, int neighbor_size) {
  if (i - neighbor_size < 0 || i + neighbor_size > rows || j - neighbor_size < 0 || j + neighbor_size > cols) {
    return 0;
  }

  float xval = 0.0f;
  float yval = 0.0f;
  for (int i1 = i - neighbor_size; i1 <= i + neighbor_size; i1++) {
    for (int i2 = j - neighbor_size; i2 <= j + neighbor_size; i2++) {
      int mask_index = (i1 - i + neighbor_size) * (2 * neighbor_size + 1) + (i2 - j + neighbor_size);
      int img_index = i1 * cols + i2;
      xval += (float) map[img_index] * sobel_x[mask_index];
      yval += (float) map[img_index] * sobel_y[mask_index];
    }
  }

  return sqrt(xval * xval + yval * yval);
}

gray *filter_mag(gray *map, int rows, int cols, int neighbor_size) {
  gray max_val = 255;
  float *gradient = malloc(rows * cols * sizeof(float));

  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      gradient[i * cols + j] = filter_pixel_mag(map, i, j, rows, cols, neighbor_size);
    }
  }

  float max_mag = 0;
  for (int i=0; i<rows * cols; i++) {
    if (gradient[i] > 0 && gradient[i] > max_mag) {
      max_mag = gradient[i];
    } else if (gradient[i] < 0 && -gradient[i] > max_mag) {
      max_mag = -gradient[i];
    }
  }

  gray *out = malloc(rows * cols * sizeof(gray));
  for (int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      out[i * cols + j] = (gray) (gradient[i * cols + j] / max_mag * max_val);
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
  int ich1, ich2, rows, cols, maxval=255, pgmraw;
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

  // gray *out = filter(graymap, rows, cols, sobel_x, 1);
  // gray *out = filter(graymap, rows, cols, sobel_y, 1);
  gray *out = filter_mag(graymap, rows, cols, 1);

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

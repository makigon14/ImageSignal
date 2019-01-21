//**********************************************************************
// Sample program: detecting a Line by Hough transform
// 
// $Id: sample.cc,v 1.1 2014/12/15 01:07:42 miura Exp miura $
//**********************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
//#include <values.h>		// for M_PI
#include <vector>
#include <complex>

using namespace std;

#define THRESHOLD (50)

//======================================================================
// Read one line from FILE **fp
//======================================================================
char *
readOneLine (char *buf, int n, FILE * fp)
{
  char *fgetsResult;

  do
    {
      fgetsResult = fgets (buf, n, fp);
    }
  while (fgetsResult != NULL && buf[0] == '#');
  // エラーや EOF ではなく、かつ、先頭が '#' の時は、次の行を読み込む

  return fgetsResult;
}

//======================================================================
// Calc angle [0, 2pi)
//======================================================================
double arg(double gx, double gy)
{
  double phi = acos(gx/sqrt(gx*gx+gy*gy));
  if(isnan(phi)) return 0.0;
  if(gy < 0) phi = 2*M_PI - phi;
  return phi;
}

//----------------------------------------------------------------------
// Image class
//----------------------------------------------------------------------
typedef int pixel_t;

class image_c
{
  // Private members
  int width;			// The number of pixels in X direction
  int height;			// The number of pixels in Y direction
  int maxValue;			// maximum value of the pixel
  pixel_t *data;		// memory area to contain the image

public:
  // Public members
  image_c(int width, int height, int maxValue);
  image_c(FILE *fp);
  inline int getWidth() const { return width; }
  inline int getHeight() const { return height; }
  inline int getMaxValue() const { return maxValue; }
  inline pixel_t &pixel(int x, int y) const { return data[x+width*y]; }
  inline pixel_t getPixel(int x, int y) const {
    if(x<0) x=0;
    if(x>=width) x=width-1;
    if(y<0) y=0;
    if(y>=height) y=height-1;
    return pixel(x, y);
  }
  inline void setPixel(int x, int y, pixel_t value) {
    if(x>=0 && x<width && y>=0 && y<height) pixel(x,y)=value;
  }
  inline pixel_t getGx(int x, int y) const {
    return -1*getPixel(x-1, y-1) +1*getPixel(x+1, y-1)
       -2*getPixel(x-1, y+0) +2*getPixel(x+1, y+0)
       -1*getPixel(x-1, y+1) +1*getPixel(x+1, y+1);
  }
  inline pixel_t getGy(int x, int y) const {
    return -1*getPixel(x-1, y-1) -2*getPixel(x+0, y-1) -1*getPixel(x+1, y-1)
      +1*getPixel(x-1, y+1) +2*getPixel(x+0, y+1) +1*getPixel(x+1, y+1);
  }
  void fill(int val);
  void drawLine(int x0, int y0, int x1, int y1, pixel_t color);
  void voteOnLine(int x0, int y0, int x1, int y1, int weight);
  void linearConv();
  void writePgmFile(FILE *fp);
  // デバッグ用出力
  void dump() const;
};

//----------------------------------------------------------------------
// image_c member functions
//----------------------------------------------------------------------

//
// image_c コンストラクタ
// サイズ指定
//
image_c::image_c(int Width, int Height, int MaxValue)
{
  width = Width;
  height = Height;
  maxValue = MaxValue;

  // メモリ領域の確保
  data = new pixel_t[width * height];

  if (data == NULL)	// メモリ確保ができなかった時はエラー
    {
      perror ("memory allocation failed");
      exit (1);
    }
}

//
// image_c constructor
// From pgm file
//
image_c::image_c (FILE * fp)
{
  char buf[128];

  // Check magic number (P5)
  if (readOneLine (buf, 128, fp) == NULL) goto error;
  if (buf[0] != 'P' || buf[1] != '5') goto error;

  // Read image size
  if (readOneLine (buf, 128, fp) == NULL) goto error;
  if (sscanf (buf, "%d %d", &width, &height) != 2) goto error;
  if (width <= 0 || height <= 0) goto error;

  // Read maximum value of the pixel
  if (readOneLine (buf, 128, fp) == NULL) goto error;
  if (sscanf (buf, "%d", &maxValue) != 1) goto error;
  if (maxValue <= 0 || maxValue >= 256) goto error;

  // Memory allocation
  data = (pixel_t *) malloc ((size_t) (sizeof(pixel_t) * width * height));

  if (data == NULL)	// Error if memory allocation failed
    {
      fputs ("out of memory\n", stderr);
      exit (1);
    }

  int i;
  for(i=0; i<width*height; i++)
    {
      int c = fgetc(fp);
      if(c == EOF) goto error;
      data[i] = c;
    }

  return;

error:
  // Error handling
  perror ("Reading PGM-RAW file was failed");
  exit (1);
}

//
// Fill the image by fixed value
//
void
image_c::fill(pixel_t val)
{
  int i;
  for(i=0; i<width*height; i++)
    data[i] = val;
}

// Draw a line
void
image_c::drawLine(int x0, int y0, int x1, int y1, pixel_t brightness)
{
  int dx = abs(x1-x0);
  int dy = abs(y1-y0);
  int sx;
  if(x0 < x1) sx = 1;
  else sx = -1;
  int sy;
  if(y0 < y1) sy = 1;
  else sy = -1;
  int err = dx-dy;
 
  while(x0!=x1 || y0!=y1)
    {
      setPixel(x0, y0, brightness);
      int e2 = 2*err;
      if(e2 > -dy)
	{
	  err = err - dy;
	  x0 = x0 + sx;
	}
      if(e2 <  dx)
	{
	  err = err + dx;
	  y0 = y0 + sy;
	}
    }
}

// Vote (draw a line)
// It is similar to drawLine() but 'weight' is added to the pixel 
void
image_c::voteOnLine(int x0, int y0, int x1, int y1, int weight)
{
  int dx = abs(x1-x0);
  int dy = abs(y1-y0);
  int sx;
  if(x0 < x1) sx = 1;
  else sx = -1;
  int sy;
  if(y0 < y1) sy = 1;
  else sy = -1;
  int err = dx-dy;
 
  while(x0!=x1 || y0!=y1)
    {
      setPixel(x0, y0, getPixel(x0, y0)+weight);
      int e2 = 2*err;
      if(e2 > -dy)
	{
	  err = err - dy;
	  x0 = x0 + sx;
	}
      if(e2 <  dx)
	{
	  err = err + dx;
	  y0 = y0 + sy;
	}
    }
}

// Linear conversion
void
image_c::linearConv()
{
  int min, max;
  min = max = pixel(0, 0);
  for(int y=0; y<height; y++) {
    for(int x=0; x<width; x++) {
      int c = pixel(x, y);
      if(c < min) min = c;
      if(max < c) max = c;
    }
  }
  for(int y=0; y<height; y++) {
    for(int x=0; x<width; x++) {
      pixel(x, y) = maxValue*(pixel(x, y)-min)/(max - min);
    }
  }
}

void
image_c::writePgmFile (FILE * fp)
{
  // Write a magic number (P5)
  if (fputs ("P5\n", fp) == EOF) goto error;

  // Write the size of the image
  if (fprintf (fp, "%d %d\n", width, height) == EOF) goto error;

  // Wreite the maximum value of the pixel
  if (fprintf (fp, "%d\n", maxValue) == EOF) goto error;

  // Write bitmap data
  int i;
  for(i=0; i<width*height; i++)
    {
      if(fputc(data[i], fp)==EOF) goto error;
    }

  return;

error:
  perror ("Writing PGM-RAW file was failed");
  exit (1);
}

void
image_c::dump() const
{
  int x, y;
  printf("   |");
  for(x=0; x<this->width; x++)
    printf("%3d", x);
  printf("\n---+");
  for(x=0; x<this->width; x++)
    printf("---");
  putchar('\n');
  for(y=0; y<this->height; y++)
    {
      printf("%3d|", y);
      for(x=0; x<this->width; x++)
	printf(" %02x", this->getPixel(x, y));

      putchar('\n');
    }
}


//**********************************************************************
// Hough transform
//**********************************************************************
#define min_m (-2.0)
#define max_m (2.0)
#define num_m (512)
#define min_c (-512.0)
#define max_c (511.0)
#define num_c (512)

void
houghTransform(image_c *result, image_c &orig)
{
  int x, y, width, height;
  
  width = result->getWidth();
  height = result->getHeight();

  for(y=0; y<height; y++)
    {
      for(x=0; x<width; x++)
	{
	  int gx = orig.getGx(x, y);
	  int gy = orig.getGy(x, y);
	  result->pixel(x, y) = (int)round(sqrt(gx*gx+gy*gy));
	}
    }

  result->linearConv();

  // Prameter space (m, c)
  // m: [min_m, max_m]
  // c: [min_c, max_c]
  image_c parameterSpace(num_m, num_c, 255);
  double delta_m = (double)(max_m - min_m)/num_m;
  double delta_c = (double)(max_c - min_c)/num_c;

  parameterSpace.fill(0);

  for(y=0; y<height; y++)
    {
      for(x=0; x<width; x++)
	{
	  // Draw a line c = (-x)*m + y in the parameter space
	  // Left point: (m, c) = (min_m, (-x)*min_m + y)
	  // Right point: (m, c) = (max_m, (-x)*max_m + y)
  	  // printf("(x,y)=(%d,%d), left: (m,c)=(%g,%g), (M,C)=(%d,%d), right: (m,c)=(%g,%g), (M,C)=(%d,%d)\n",
	  // 	 x, y,
	  // 	 min_m, (-x)*min_m + y,
	  // 	 0, (int)round(((-x)*min_m + y - min_c)/delta_c),
	  // 	 max_m, (-x)*max_m + y,
	  // 	 num_m-1, (int)round(((-x)*max_m + y - min_c)/delta_c));
	  if(result->pixel(x, y) > THRESHOLD) {
	    parameterSpace.voteOnLine( 0, (int)round(((-x)*min_m + y - min_c)/delta_c),
				       num_m, (int)round(((-x)*max_m + y - min_c)/delta_c),
				       orig.getPixel(x, y) );
	  }
	}
    }

  // parameterSpace.linearConv();
  // parameterSpace.writePgmFile(stdout);
  // exit(0);
	
  int max = 0;
  int max_M;
  int max_C;
  for(int C=0; C<parameterSpace.getHeight(); C++)
    {
      for(int M=0; M<parameterSpace.getWidth(); M++)
	{
	  if(max < parameterSpace.getPixel(M, C))
	    {
	      max = parameterSpace.getPixel(M, C);
	      max_C = C;
	      max_M = M;
	    }
	}
    }

  // copy image
  for(y=0; y<height; y++)
    {
      for(x=0; x<width; x++)
	{
	  result->setPixel(x, y, orig.getPixel(x, y));
	}
    }

  // Draw the detected line (y = m * x + c) in black
  // Left point: (x, y) = (0, c)
  // Right point: (x, y) = (width-1, m*(width-1)+c)
  double m = max_M * delta_m + min_m;
  double c = max_C * delta_c + min_c;
  fprintf(stderr, "(m, c)=(%g, %g)\n", m, c);
  result->drawLine(0, (int)round(c),
		   width-1, (int)round(m*(width-1)+c),
		   0);

  return;
}

//**********************************************************************
// メイン
//**********************************************************************
int
main (int argc, char **argv)
{
  image_c orig (stdin);

  image_c result (orig.getWidth(), orig.getHeight(), orig.getMaxValue());

  houghTransform(&result, orig);

  result.writePgmFile(stdout);

  return 0;
}

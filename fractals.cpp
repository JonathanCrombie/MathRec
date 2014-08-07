/* Public Domain.  See the LICENSE file. */

/* This is a command line program to generate the    */
/* Mandelbrot Set and Julia Sets.                    */
 
/* There are no dependencies.  To compile on linux,  */
/* try: gcc fractals.cpp -o fractals                 */

/* Using Visual C++ in Windows, the following        */
/* worked from a command prompt: cl fractals.cpp     */

#include <stdio.h>
#include <stdlib.h>

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <fcntl.h>
#include <io.h>
#endif

#include <math.h>
#include <string.h>
#include <ctype.h>

struct pixel
{
    unsigned char   red;
    unsigned char   green;
    unsigned char   blue;
};

void printusage();
int Get2Tuple( char*, double*, double* );
int Get2Tuple( char*, long*, long* );
void initpal(struct pixel *);

const char* VersionStr = "1.0.1";
const unsigned char CRLF[2] = {0x0D,0x0A};

int main( int argc, char* argv[] ) {

#if defined(_WIN32) && !defined(__CYGWIN__)
 if ( _setmode( _fileno( stdout ), _O_BINARY ) == -1 ) {
    printf( "Cannot set stdout to binary mode.  Exiting." );
    return -1;
 }
#endif

  char*     userfilename = NULL;
  int       user_capk    = -1;
  double    user_centerx = 0.0;
  double    user_centery = 0.0;
  int       user_centeroverride = 0;
  double    user_julia_r = 0.0;
  double    user_julia_i = 0.0;
  int       MakeJuliaSet = 0;
  long      user_resolx = 0.0;
  long      user_resoly = 0.0;
  int       user_resolutionoverride = 0;
  double    user_zoomlevel = -1.0;

  long i;
  for ( i = 1; i < argc; ) {
    int ishyphen = argv[i][0] == '-';
    long len = strlen( argv[i] );
    long nextlen = (i+1) < argc ? strlen( argv[i+1] ) : -1;

    int argsprocessed = 0;
    if ( ishyphen && len > 1 )
      argsprocessed = 1;

    if ( argsprocessed >= 1 ) {
      char useroption = argv[i][1];
      char* optionvalue = NULL;
      if ( len >= 3 )
        optionvalue = &argv[i][2];

      switch ( useroption ) {
       case 'c':  // center point  (x,y)
        if ( optionvalue == NULL && nextlen > 0 ) {
          optionvalue = argv[i+1];
          argsprocessed = 2;
        }
        if ( optionvalue != NULL )
          user_centeroverride = !Get2Tuple( optionvalue, &user_centerx, &user_centery );
        break;
       case 'h':
        printusage();
        return 0;
        break;
       case 'j':  // julia set constant value (the real part and imaginary part)
        if ( optionvalue == NULL && nextlen > 0 ) {
          optionvalue = argv[i+1];
          argsprocessed = 2;
        }
        if ( optionvalue != NULL )
          MakeJuliaSet = !Get2Tuple( optionvalue, &user_julia_r, &user_julia_i );
        break;
       case 'm':  // maximum number of iterations per pixel
        if ( optionvalue == NULL && nextlen > 0 ) {
          optionvalue = argv[i+1];
          argsprocessed = 2;
        }
        if ( optionvalue != NULL )
          user_capk = abs(atoi( optionvalue ));
        break;
       case 'o':  // output file name
        if ( optionvalue == NULL && nextlen > 0 ) {
          optionvalue = argv[i+1];
          argsprocessed = 2;
        }
        if ( optionvalue != NULL )
          userfilename = strdup( optionvalue );
        break;
       case 'r':  // image resolution
        if ( optionvalue == NULL && nextlen > 0 ) {
          optionvalue = argv[i+1];
          argsprocessed = 2;
        }
        if ( optionvalue != NULL )
          user_resolutionoverride = !Get2Tuple( optionvalue, &user_resolx, &user_resoly );
        break;
       case 'v':
        printf( "fractals version %s\n", VersionStr );
        return 0;
        break;
       case 'z':
        if ( optionvalue == NULL && nextlen > 0 ) {
          optionvalue = argv[i+1];
          argsprocessed = 2;
        }
        if ( optionvalue != NULL )
          user_zoomlevel = fabs(atof( optionvalue ));
        break;
       default:
        break;
      }
    }

    if ( argsprocessed == 2 )
      i += 2;
    else
      i++;
  }

  double centerx = 0.0;
  double centery = 0.0;
  if ( !MakeJuliaSet )  // ie. Make the Mandelbrot Set
    centerx = -0.75;
  if ( user_centeroverride ) {
    centerx = user_centerx;
    centery = user_centery;
  }

  double c_r = 0.0;
  double c_i = 0.0;
  if ( MakeJuliaSet ) {
    c_r = user_julia_r;
    c_i = user_julia_i;
  }

  int capk = 2048;
  if ( user_capk > 0 && user_capk < 10000000 )
    capk = user_capk;

  FILE* fpout = stdout;
  if ( userfilename != NULL ) {
    FILE* fdtest = fopen( userfilename, "r" );
    if ( fdtest != NULL ) {
      printf("Output file \"%s\" already exists.  Refusing to overwrite.  Exiting.\n\n", userfilename );
      fclose( fdtest );
      free( userfilename );
      return -1;
    }
    fpout = fopen( userfilename, "wb" );
    if ( fpout == NULL ) {
      printf("Error: Could not open file \"%s\" for write.  Exiting.\n\n", userfilename );
      free( userfilename );
      return -1;
    }
  }

  long resolx = 1024;  // horizontal resolution in pixels
  long resoly = 768;   // vertical resolution in pixels
  if ( user_resolutionoverride ) {
    resolx = user_resolx;
    resoly = user_resoly;
  }

  double zoomlevel = 1.0;  // zoomlevel of 1.0 arbitrarily defined to be an x-width of 3.1.
  if ( user_zoomlevel > 0.00001 && user_zoomlevel < 10000000.0 )
    zoomlevel = user_zoomlevel;
  double fulldx = 3.1 / zoomlevel;
  double fulldy = (3.1 / zoomlevel) * ((double)resoly /(double)resolx);

  double pixelwidth = fulldx/(double)resolx;  // Could of just as easily been fulldy/(double)resoly
  double halfpixel = pixelwidth / 2.0;

  double xmin = centerx - fulldx / 2.0;
  double xminplushalf = xmin + halfpixel; // like to use the middles of pixels
  double ymax = centery + fulldy / 2.0;
  double ymaxlesshalf = ymax - halfpixel; // like to use the middles of pixels


  const double m  = 100.0;  // min norm to be considered an escapee


  struct pixel holdpal[256];
  initpal( holdpal );

  fprintf( fpout, "P6" );
  fwrite( CRLF, 1, 2, fpout );

  fprintf( fpout, "%ld %ld", resolx, resoly );
  fwrite( CRLF, 1, 2, fpout );

  fprintf( fpout, "255" );
  fwrite( CRLF, 1, 2, fpout );

  long x,y;
  for ( y = 0; y < resoly; y++ ) {
    for ( x = 0; x < resolx; x++ ) {
      double z_r = 0.0;
      double z_i = 0.0;

      if ( MakeJuliaSet ) {
        z_r = xminplushalf + x * pixelwidth;
        z_i = ymaxlesshalf - y * pixelwidth;
      }
      else {  // Make the Mandelbrot Set
        c_r = xminplushalf + x * pixelwidth;
        c_i = ymaxlesshalf - y * pixelwidth;
      }

      int k = -1;
      double norm = 0.0;

      double z_r_save = z_r;
      while ( norm < m && k < capk ) {  // repeatedly iterating z = z^2 + c  where z & c are complex numbers
        z_r_save = z_r;
        z_r = z_r_save * z_r_save - z_i * z_i + c_r;
        z_i = 2 * z_r_save * z_i + c_i;
        k++;
        norm = z_r * z_r + z_i * z_i;
      }

      if ( k == capk )
          k = 255;
      else
          k %= 254;

      fwrite( &holdpal[k], 1, 3, fpout );
    }
  }

  if ( fpout != stdout ) {
    fclose(fpout);
    fpout = NULL;
  }

  if ( userfilename != NULL ) {
    free( userfilename );
    userfilename = NULL;
  }

  return 0;
}

void printusage() {
  printf( "\n" );
  printf( "fractals version %s\n\n", VersionStr );

  printf( "usage: fractals [options]\n\n" );

  printf( "options:\n" );
  printf( "  -c real_x,real_y    -- specifies the center coordinates (real_x,real_y).\n" );
  printf( "  -h                  -- prints this help and exits.\n" );
  printf( "  -j p,q              -- generate a Julia Set with complex c = p + qi.\n" );
  printf( "  -m integer          -- specifies the maximum # of iterations per pixel\n");
  printf( "                         before stopping.\n" );
  printf( "  -o filename         -- save to this output file.\n" );
  printf( "  -r integer,integer  -- image resolution.\n" );
  printf( "  -v                  -- print version and exit.\n" );
  printf( "  -z real             -- set zoom level to real.\n\n" );

  printf( " modes:\n" );
  printf( "   fractals has 2 modes.  The Mandelbrot mode is the default, but it will\n" );
  printf( "   switch to Julia Set mode if a \"-j p,q\" option is used.\n\n" );

  printf( " defaults:\n" );
  printf( "   -- The default center is (0.75,0.0) for Mandelbrot mode and (0.0,0.0) for\n" );
  printf( "      Julia Set mode.\n" );
  printf( "   -- The default for m is 2048.\n" );
  printf( "   -- The default output is to stdout.\n" );
  printf( "   -- The default image resolution is 1024x768.\n" );
  printf( "   -- The default zoom level is 1.0 which is a real x-width of 3.1.\n\n" );

  printf( " examples:\n" );
  printf( "   fractals > mset.ppm\n" );
  printf( "     -- produces a Mandelbrot Set called \"mset.ppm\".\n" );
  printf( "   fractals -o mset.ppm\n" );
  printf( "     -- same result as \"fractals > mset.ppm\".\n" );
  printf( "   fractals | pnmtopng > mset.png\n" );
  printf( "     -- create a loss-less compressed .png file \"mset.png\".  Need \"netpbm\"\n" );
  printf( "        installed.\n" );
  printf( "   fractals | pnmtojpeg > mset.jpg\n" );
  printf( "     -- create a lossy compressed jpeg file \"mset.jpg\".  Need \"netpbm\"\n" );
  printf( "        installed.\n" );
  printf( "   fractals -j -.194,.6557 > jset.ppm\n" );
  printf( "     -- create the Julia Set with c = -.194 + .6557i and save in \"jset.ppm\".\n" );
  printf( "   fractals -j-.194,.6557 -c-.32,0.27 -r1280x960 -m3000 -z4.75 > jset2.ppm\n" );
  printf( "     -- create the Julia Set with c = -.194 + .6557i and save in \"jset2.ppm\".\n" );
  printf( "        set center to (-0.32,0.27), resolution to 1280 by 960 pixels, max\n" );
  printf( "        iterations to 3000 and zoom level to 4.75.\n\n" );

  printf( "\n\n" );
}

// parse out two doubles from inputstr
int Get2Tuple( char* inputstr, double* first, double* second ) {

  char* tempstr = strdup( inputstr );

  int i;
  int len = strlen( tempstr );

  char* begdouble1 = NULL;
  char* begdouble2 = NULL;

  int findbegin = 1;
  int findnumberparts = 0;
  for ( i = 0; i < len; i++ ) {
    int isnumberpart = isdigit( tempstr[i] ) ||  tempstr[i] == '-' || tempstr[i] == '.';
    if ( findbegin && isnumberpart ) {
      if ( begdouble1 == NULL )
        begdouble1 = tempstr + i;
      else
        begdouble2 = tempstr + i;
      findbegin = 0;
      findnumberparts = 1;
    }
    else if ( findnumberparts && !isnumberpart ) {
      findnumberparts = 0;
      tempstr[i] = '\0';
      if ( begdouble2 == NULL )
        findbegin = 1;
    }
  }

  int fail = 1;
  if ( begdouble1 != NULL && begdouble2 != NULL ) {
    *first  = atof( begdouble1 );
    *second = atof( begdouble2 );
    fail = 0;
  }

  free( tempstr );
  tempstr = NULL;

  return fail;
}

// parse out two longs from inputstr
int Get2Tuple( char* inputstr, long* first, long* second ) {

  char* tempstr = strdup( inputstr );

  int i;
  int len = strlen( tempstr );

  char* beglong1 = NULL;
  char* beglong2 = NULL;

  int findbegin = 1;
  int findnumberparts = 0;
  for ( i = 0; i < len; i++ ) {
    int isnumberpart = isdigit( tempstr[i] );
    if ( findbegin && isnumberpart ) {
      if ( beglong1 == NULL )
        beglong1 = tempstr + i;
      else
        beglong2 = tempstr + i;
      findbegin = 0;
      findnumberparts = 1;
    }
    else if ( findnumberparts && !isnumberpart ) {
      findnumberparts = 0;
      tempstr[i] = '\0';
      if ( beglong2 == NULL )
        findbegin = 1;
    }
  }

  int fail = 1;
  if ( beglong1 != NULL && beglong2 != NULL ) {
    *first  = atol( beglong1 );
    *second = atol( beglong2 );
    fail = 0;
  }

  free( tempstr );
  tempstr = NULL;

  return fail;
}

// create a palette
void initpal( struct pixel holdpal[256] ) {
  int         i;

  for (i = 0; i < 64; i++) {
    holdpal[i].red = 125 - i;
    holdpal[i].green = 61 + i;
    holdpal[i].blue = 254 - (i * 2);
  }

  for (i = 64; i < 128; i++) {
    holdpal[i].red = 61 + (i - 64);
    holdpal[i].green = 125 + ((i - 64) * 2);
    holdpal[i].blue = 125 - (i - 64);
  }

  for (i = 128; i < 192; i++) {
    holdpal[i].red = 125 + ((i - 128) * 2);
    holdpal[i].green = 254 - ((i - 128) * 2);
    holdpal[i].blue = 61 + (i - 128);
  }

  for (i = 192; i < 255; i++) {
    holdpal[i].red = 254 - ((i - 192) * 2);
    holdpal[i].green = 125 - (i - 192);
    holdpal[i].blue = 125 + ((i - 192) * 2);
  }

  holdpal[255].red = 0;
  holdpal[255].green = 0;
  holdpal[255].blue = 0;
}


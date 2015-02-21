#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "FreeImage.h"

float matrix[4][4][3] ;             /* Matrix for quadratic form */

float sinc(float x) {               /* Supporting sinc function */
  if (fabs(x) < 1.0e-4) return 1.0 ;
  else return(sin(x)/x) ;
}

struct squareFloatImgBlob {
    int width;
    float* data;
};
typedef struct squareFloatImgBlob squareFloatImgBlob;


void loadFromFloatFormat(const char * filename, squareFloatImgBlob* outBlob) {
  FILE *fp ;
  assert(fp = fopen(filename,"rb")) ;

  fseek(fp, 0, SEEK_END);

  long fileSize = ftell(fp);

  int width = sqrt(fileSize/12);

  fseek(fp, 0, SEEK_SET);

  outBlob->width = width;
  outBlob->data = (float*) malloc(fileSize);

  size_t bytesRead = fread(outBlob->data, 1, fileSize, fp);

  printf("File %s has width %d and %d\n", filename, width, bytesRead);

  fclose(fp) ;
}

int min(int a, int b) {
    return (a < b ? a : b);
}

void squareFloatImgBlob_init(squareFloatImgBlob* blob, int width) {
    blob->width = width;
    blob->data = malloc(width*width*3*sizeof(float));
}

void prefilterAngleMap(const squareFloatImgBlob* inFloatBlob, float coeffs[9][3]) {

  /* The main integration routine.  Of course, there are better ways
     to do quadrature but this suffices.  Calls updatecoeffs to
     actually increment the integral. Width is the size of the
     environment map */

    int width = inFloatBlob->width;

    float* hdr = inFloatBlob->data;


    for (int i = 0 ; i < width ; i++) {
        for (int j = 0 ; j < width ; j++) {

            /* We now find the cartesian components for the point (i,j) */
            float u,v,r,theta,phi,x,y,z,domega ;

            v = (width/2.0 - i)/(width/2.0);  /* v ranges from -1 to 1 */
            u = (j - width/2.0)/(width/2.0);    /* u ranges from -1 to 1 */
            r = sqrtf(u*u + v*v) ;               /* The "radius" */

            if (r <= 1.0) {           /* Consider only circle with r<1 */

                theta = M_PI*r ;                    /* theta parameter of (i,j) */
                phi = atan2(v,u) ;                /* phi parameter */

                x = sin(theta)*cos(phi) ;         /* Cartesian components */
                y = sin(theta)*sin(phi) ;
                z = cos(theta) ;

                domega = (2*M_PI/width)*(2*M_PI/width)*sinc(theta) ;

                for (int channel = 0 ; channel < 3 ; channel++) {
                    float c ; /* A different constant for each coefficient */

                    /* L_{00}.  Note that Y_{00} = 0.282095 */
                    c = 0.282095 ;
                    coeffs[0][channel] += hdr[i * width * 3 + j * 3 + channel]*c*domega ;

                    /* L_{1m}. -1 <= m <= 1.  The linear terms */
                    c = 0.488603 ;
                    coeffs[1][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*y)*domega ;   /* Y_{1-1} = 0.488603 y  */
                    coeffs[2][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*z)*domega ;   /* Y_{10}  = 0.488603 z  */
                    coeffs[3][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*x)*domega ;   /* Y_{11}  = 0.488603 x  */

                    /* The Quadratic terms, L_{2m} -2 <= m <= 2 */

                    /* First, L_{2-2}, L_{2-1}, L_{21} corresponding to xy,yz,xz */
                    c = 1.092548 ;
                    coeffs[4][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*x*y)*domega ; /* Y_{2-2} = 1.092548 xy */
                    coeffs[5][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*y*z)*domega ; /* Y_{2-1} = 1.092548 yz */
                    coeffs[7][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*x*z)*domega ; /* Y_{21}  = 1.092548 xz */

                    /* L_{20}.  Note that Y_{20} = 0.315392 (3z^2 - 1) */
                    c = 0.315392 ;
                    coeffs[6][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*(3*z*z-1))*domega ;

                    /* L_{22}.  Note that Y_{22} = 0.546274 (x^2 - y^2) */
                    c = 0.546274 ;
                    coeffs[8][channel] += hdr[i * width * 3 + j * 3 + channel]*(c*(x*x-y*y))*domega ;
                }
            }
        }

    }
}

void tomatrix(float coeffs[9][3]) {

  /* Form the quadratic form matrix (see equations 11 and 12 in paper) */

  int col ;
  float c1,c2,c3,c4,c5 ;
  c1 = 0.429043 ; c2 = 0.511664 ;
  c3 = 0.743125 ; c4 = 0.886227 ; c5 = 0.247708 ;

  for (col = 0 ; col < 3 ; col++) { /* Equation 12 */

    matrix[0][0][col] = c1*coeffs[8][col] ; /* c1 L_{22}  */
    matrix[0][1][col] = c1*coeffs[4][col] ; /* c1 L_{2-2} */
    matrix[0][2][col] = c1*coeffs[7][col] ; /* c1 L_{21}  */
    matrix[0][3][col] = c2*coeffs[3][col] ; /* c2 L_{11}  */

    matrix[1][0][col] = c1*coeffs[4][col] ; /* c1 L_{2-2} */
    matrix[1][1][col] = -c1*coeffs[8][col]; /*-c1 L_{22}  */
    matrix[1][2][col] = c1*coeffs[5][col] ; /* c1 L_{2-1} */
    matrix[1][3][col] = c2*coeffs[1][col] ; /* c2 L_{1-1} */

    matrix[2][0][col] = c1*coeffs[7][col] ; /* c1 L_{21}  */
    matrix[2][1][col] = c1*coeffs[5][col] ; /* c1 L_{2-1} */
    matrix[2][2][col] = c3*coeffs[6][col] ; /* c3 L_{20}  */
    matrix[2][3][col] = c2*coeffs[2][col] ; /* c2 L_{10}  */

    matrix[3][0][col] = c2*coeffs[3][col] ; /* c2 L_{11}  */
    matrix[3][1][col] = c2*coeffs[1][col] ; /* c2 L_{1-1} */
    matrix[3][2][col] = c2*coeffs[2][col] ; /* c2 L_{10}  */
    matrix[3][3][col] = c4*coeffs[0][col] - c5*coeffs[6][col] ;
                                            /* c4 L_{00} - c5 L_{20} */
  }
}

void printMatricesToGlslDeclaration() {
    printf("\nconst mat4 gracered = mat4(\n") ;
    for (int i = 0 ; i < 4 ; i++) {
        printf("    ");
        for (int j = 0 ; j < 4 ; j++) {
            if (i == 3 && j == 3) {
                printf("%9.6f ", matrix[i][j][0]) ;
            }
            else {
                printf("%9.6f, ", matrix[i][j][0]) ;
            }
        }
        printf("\n") ;
    }
    printf(") ;\n");

    printf("\nconst mat4 gracegreen = mat4(\n") ;
    for (int i = 0 ; i < 4 ; i++) {
        printf("    ");
        for (int j = 0 ; j < 4 ; j++) {
            if (i == 3 && j == 3) {
                printf("%9.6f ", matrix[i][j][1]) ;
            }
            else {
                printf("%9.6f, ", matrix[i][j][1]) ;
            }
        }
        printf("\n") ;
    }
    printf(") ;\n");

    printf("\nconst mat4 graceblue = mat4(\n") ;
    for (int i = 0 ; i < 4 ; i++) {
        printf("    ");
        for (int j = 0 ; j < 4 ; j++) {
            if (i == 3 && j == 3) {
                printf("%9.6f ", matrix[i][j][2]) ;
            }
            else {
                printf("%9.6f, ", matrix[i][j][2]) ;
            }
        }
        printf("\n") ;
    }
    printf(") ;\n");
}

void swapFloat(float* lhs, float* rhs) {
    float tmp = *lhs;

    *lhs = *rhs;
    *rhs = tmp;
}


//void squareByteImgBlob_loadFromDib(squareByteImgBlob* blob, FIBITMAP* dib) {
//    squareByteImgBlob_init(blob, FreeImage_GetWidth(dib));

//    FreeImage_ConvertToRawBits(blob->data, dib, FreeImage_GetWidth(dib) * 3, 24, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);

//    // Swapping blue and red channel since FIF is BGR and masks seems not to work
//    for (int i = 0; i < blob->width; ++i) {
//        for (int j = 0; j < blob->width; ++j) {
//            swapByte(&(blob->data[i * blob->width * 3 + j * 3 + 0]), &(blob->data[i * blob->width * 3 + j * 3 + 2]));
//        }
//    }
//}

char* getFilenameWithoutExt(const char* filename) {
    const char *dot = strrchr(filename, '.');

    if(!dot || dot == filename) return "";

    size_t size = dot - filename;

    char* ret = malloc(size + 1);

    memcpy(ret, filename, size);
    ret[size] = 0;

    return ret;
}

char* getFilenameExt(const char* filename) {
    const char *dot = strrchr(filename, '.');

    if(!dot || dot == filename) return "";

    return strdup(dot + 1);
}


int main( int argc, char *argv[] )
{
    if(argc < 2)
    {
        printf("Il n'y pas de fichier en entrée.\n");
        exit(-1);
    }

    const char* file = argv[1];
    char * ext;
    if( ( ext = strrchr( file, '.') ) != NULL )
    {
        printf("Usage : Ne pas renseigner l'extension de la probe.\n");
        exit(-2);
    }

    char file_float[ strlen(file) + 6 ];
    strcpy( file_float, file ); 
    strcat( file_float, ".float" ); 

    char filename[ strlen(file) + 4 ];
    strcpy( filename, file ); 
    strcat( filename, ".hdr" ); 
    
    squareFloatImgBlob floatProbe;
    
    loadFromFloatFormat( file_float, &floatProbe);
    
    //const char* filename = "grace_probe.hdr";

    /*
    FI_ENUM(FREE_IMAGE_FORMAT) {
        FIF_UNKNOWN = -1,
        FIF_BMP		= 0,
        FIF_ICO		= 1,
        FIF_JPEG	= 2,
        FIF_JNG		= 3,
        FIF_KOALA	= 4,
        FIF_LBM		= 5,
        FIF_IFF = FIF_LBM,
        FIF_MNG		= 6,
        FIF_PBM		= 7,
        FIF_PBMRAW	= 8,
        FIF_PCD		= 9,
        FIF_PCX		= 10,
        FIF_PGM		= 11,
        FIF_PGMRAW	= 12,
        FIF_PNG		= 13,
        FIF_PPM		= 14,
        FIF_PPMRAW	= 15,
        FIF_RAS		= 16,
        FIF_TARGA	= 17,
        FIF_TIFF	= 18,
        FIF_WBMP	= 19,
        FIF_PSD		= 20,
        FIF_CUT		= 21,
        FIF_XBM		= 22,
        FIF_XPM		= 23,
        FIF_DDS		= 24,
        FIF_GIF     = 25,
        FIF_HDR		= 26,
        FIF_FAXG3	= 27,
        FIF_SGI		= 28,
        FIF_EXR		= 29,
        FIF_J2K		= 30,
        FIF_JP2		= 31,
        FIF_PFM		= 32,
        FIF_PICT	= 33,
        FIF_RAW		= 34,
        FIF_WEBP	= 35,
        FIF_JXR		= 36
    };
    */

    FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
    fif = FreeImage_GetFileType(filename, 0);

    printf("Format is %d \n", fif);

    if(fif == FIF_UNKNOWN)
    {
       printf("Cannot get filetype from signature, trying extension\n");
       // If we can't get the signature, try to guess the file format from the file extension
       fif = FreeImage_GetFIFFromFilename(filename);
    }

    if(fif == FIF_UNKNOWN)
    {
       printf("Cannot get filetype from extension neither, aborting \n");

       exit(EXIT_FAILURE);
    }

    FIBITMAP* dib = FreeImage_Load(fif, filename, 0);

//    FI_ENUM(FREE_IMAGE_TYPE) {
//        FIT_UNKNOWN = 0,	// unknown type
//        FIT_BITMAP  = 1,	// standard image			: 1-, 4-, 8-, 16-, 24-, 32-bit
//        FIT_UINT16	= 2,	// array of unsigned short	: unsigned 16-bit
//        FIT_INT16	= 3,	// array of short			: signed 16-bit
//        FIT_UINT32	= 4,	// array of unsigned long	: unsigned 32-bit
//        FIT_INT32	= 5,	// array of long			: signed 32-bit
//        FIT_FLOAT	= 6,	// array of float			: 32-bit IEEE floating point
//        FIT_DOUBLE	= 7,	// array of double			: 64-bit IEEE floating point
//        FIT_COMPLEX	= 8,	// array of FICOMPLEX		: 2 x 64-bit IEEE floating point
//        FIT_RGB16	= 9,	// 48-bit RGB image			: 3 x 16-bit
//        FIT_RGBA16	= 10,	// 64-bit RGBA image		: 4 x 16-bit
//        FIT_RGBF	= 11,	// 96-bit RGB float image	: 3 x 32-bit IEEE floating point
//        FIT_RGBAF	= 12	// 128-bit RGBA float image	: 4 x 32-bit IEEE floating point
//    };

    FREE_IMAGE_TYPE type = FreeImage_GetImageType(dib);

    printf("Type is %d \n", type);

    // Allocate a raw buffer
    int width = FreeImage_GetWidth(dib);
    int height = FreeImage_GetHeight(dib);
    int pitch = FreeImage_GetPitch(dib);
    int bpp = FreeImage_GetBPP(dib);
//    BYTE *bits = (BYTE*)malloc(height * pitch);

    squareFloatImgBlob floatProbe2;

//    loadFromFloatFormat("grace_probe.float", &floatProbe2);

    squareFloatImgBlob_init(&floatProbe2, width);

    // convert the bitmap to raw bits (top-left pixel first)
//    FreeImage_ConvertToRawBits((BYTE*) floatProbe2.data, dib, FreeImage_GetWidth(dib) * 12, bpp, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);
//    FreeImage_Unload(dib);

//    for (int i = 0; i < floatProbe2.width; ++i) {
//        for (int j = 0; j < floatProbe2.width; ++j) {
//            swapFloat(&(floatProbe2.data[i * floatProbe2.width * 3 + j * 3 + 0]), &(floatProbe2.data[i * floatProbe2.width * 3 + j * 3 + 2]));
//        }
//    }

//    floatProbe2.data = (float*) FreeImage_GetBits(dib);

    FIRGBF testF;

    float coeffs[9][3] = {0} ;

    for(int y = 0; y < FreeImage_GetHeight(dib); y++) {
        FIRGBF *bits = (FIRGBF *)FreeImage_GetScanLine(dib, y);

        for(int x = 0; x < FreeImage_GetWidth(dib); x++) {

            int i = y;
            int j = x;

            floatProbe2.data[i * floatProbe2.width * 3 + j * 3 + 0] = bits[x].red;
            floatProbe2.data[i * floatProbe2.width * 3 + j * 3 + 1] = bits[x].green;
            floatProbe2.data[i * floatProbe2.width * 3 + j * 3 + 2] = bits[x].blue;
        }
    }


    for (int i = 500; i < 501; ++i) {
        for (int j = 500; j < 505; ++j) {
            printf("%f  ", floatProbe.data[i * floatProbe.width * 3 + j * 3 + 1]);
        }

        printf("\n");
    }

    for (int i = 500; i < 501; ++i) {
        for (int j = 500; j < 505; ++j) {
            printf("%f  ", floatProbe2.data[i * floatProbe2.width * 3 + j * 3 + 1]);
        }

        printf("\n");
    }

    FIRGBF* sl = (FIRGBF*) FreeImage_GetScanLine(dib, 500);

    printf("From scanline : %f\n", sl[500].green);

    printf("%f  \n", floatProbe.data[499 * floatProbe.width * 3 + 499 * 3 + 1]);
    printf("%f  \n", floatProbe2.data[500 * floatProbe2.width * 3 + 500 * 3 + 1]);

    prefilterAngleMap(&floatProbe, coeffs);

    printf("Filename            : %s\n", filename);
    printf("Filename (no ext)   : %s\n", getFilenameWithoutExt(filename));
    printf("Extension           : %s\n", getFilenameExt(filename));

    char* binaryFilename = malloc(strlen(filename) + 1 + strlen("leadrshc") + 1);
    strcpy(binaryFilename, getFilenameWithoutExt(filename));
    strcat(binaryFilename, ".leadrshc");

    printf("\n%s\n", binaryFilename);

    FILE* fp = NULL;

    fp = fopen(binaryFilename, "wb");

    fwrite(coeffs, sizeof(float), 9*3, fp);

    fclose(fp);

    for (int i = 0; i < 9*3; ++i) {
        (coeffs[0])[i] = 0.f;
    }

    fp = NULL;

    fp = fopen(binaryFilename, "rb");

    fread(coeffs, sizeof(float), 9*3, fp);

    fclose(fp);


//    printf("\n         Lighting Coefficients\n\n") ;
//    printf("(l,m)       RED        GREEN     BLUE\n") ;

//    printf("L_{0,0}   %9.6f %9.6f %9.6f\n",
//       coeffs[0][0],coeffs[0][1],coeffs[0][2]) ;
//    printf("L_{1,-1}  %9.6f %9.6f %9.6f\n",
//       coeffs[1][0],coeffs[1][1],coeffs[1][2]) ;
//    printf("L_{1,0}   %9.6f %9.6f %9.6f\n",
//       coeffs[2][0],coeffs[2][1],coeffs[2][2]) ;
//    printf("L_{1,1}   %9.6f %9.6f %9.6f\n",
//       coeffs[3][0],coeffs[3][1],coeffs[3][2]) ;
//    printf("L_{2,-2}  %9.6f %9.6f %9.6f\n",
//       coeffs[4][0],coeffs[4][1],coeffs[4][2]) ;
//    printf("L_{2,-1}  %9.6f %9.6f %9.6f\n",
//       coeffs[5][0],coeffs[5][1],coeffs[5][2]) ;
//    printf("L_{2,0}   %9.6f %9.6f %9.6f\n",
//       coeffs[6][0],coeffs[6][1],coeffs[6][2]) ;
//    printf("L_{2,1}   %9.6f %9.6f %9.6f\n",
//       coeffs[7][0],coeffs[7][1],coeffs[7][2]) ;
//    printf("L_{2,2}   %9.6f %9.6f %9.6f\n",
//       coeffs[8][0],coeffs[8][1],coeffs[8][2]) ;

//    tomatrix(coeffs);

//    printMatricesToGlslDeclaration();


    return 0;
}


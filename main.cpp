#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>

#include "IL/il.h"

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

float floatReverseEndian(float val) {
    unsigned char* pSrc = (unsigned char*) &val;

    float dst = 0;
    unsigned char* pDst = (unsigned char*) &dst;

    for (int k = 0; k < 4; ++k) {
        pDst[k] = pSrc[3 - k];
    }

    return dst;
}

void loadFromFloatFormat(const char * filename, squareFloatImgBlob* outBlob) {
  FILE *fp ;
  assert(fp = fopen(filename,"rb")) ;

  fseek(fp, 0, SEEK_END);

  long fileSize = ftell(fp);

  int width = sqrt(fileSize/12);

  fseek(fp, 0, SEEK_SET);

  outBlob->width = width;
  outBlob->data = (float*) malloc(fileSize);

//  for (int i = 0; i < width*width*3; ++i) {
//      float val;

//      fread(&val, sizeof(float), 1, fp);

//      outBlob->data[i] = floatReverseEndian(val);
//  }

  size_t bytesRead = fread(outBlob->data, 1, fileSize, fp);

//  printf("File %s has width %d and %d\n", filename, width, bytesRead);

  fclose(fp) ;
}

void squareFloatImgBlob_init(squareFloatImgBlob* blob, int width) {
    blob->width = width;
    blob->data = (float*) malloc(width*width*3*sizeof(float));
}

void squareFloatImgBlob_reverseEndian(squareFloatImgBlob* blob) {
//    for (int i = 0; i < blob->width*blob->width*3; ++i) {
//        blob->data[i] = floatReverseEndian(blob->data[i]);
//    }
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


char* getFilenameWithoutExt(const char* filename) {
    const char *dot = strrchr(filename, '.');

    if(!dot || dot == filename) return (char*) "";

    size_t size = dot - filename;

    char* ret = (char*) malloc(size + 1);

    memcpy(ret, filename, size);
    ret[size] = 0;

    return ret;
}

char* getFilenameExt(const char* filename) {
    const char *dot = strrchr(filename, '.');

    if(!dot || dot == filename) return (char*) "";

    return strdup(dot + 1);
}


int main(int argc, char** argv)
{
    squareFloatImgBlob floatProbe;

    ilInit();

    if(argc != 2)
    {
        std::cerr << "Wrong arguments number.\n";
        std::cerr << "Usage : IrradianceMap pic_filename \n";

        exit(EXIT_FAILURE);
    }

    const char* filename = argv[1];

    if (strcmp(getFilenameExt(filename), "float") == 0) {
        loadFromFloatFormat(filename, &floatProbe);
    }
    else {
        unsigned int width;
        unsigned int height;
        float* data;
        unsigned int bytesPerPixel;
        int format;
        int type;
        int numChannels;

        ILuint imageId;

        bool reversed = false;

        // The image name to return.
        ilGenImages(1, &imageId); // Grab a new image name.
        ilBindImage(imageId);
        ilEnable(IL_ORIGIN_SET);

        if (reversed)
            ilOriginFunc(IL_ORIGIN_UPPER_LEFT);
        else
            ilOriginFunc(IL_ORIGIN_LOWER_LEFT);

        if (ilLoadImage(filename) != IL_TRUE) {
            std::cerr << "Error loading image from: " << filename << '\n';
            return false;
        }

        data = (float*) ilGetData();

        width = ilGetInteger(IL_IMAGE_WIDTH);
        height = ilGetInteger(IL_IMAGE_HEIGHT);

        bytesPerPixel = ilGetInteger(IL_IMAGE_BYTES_PER_PIXEL);

        format = ilGetInteger(IL_IMAGE_FORMAT);
        type = ilGetInteger(IL_IMAGE_TYPE);
        numChannels = ilGetInteger(IL_IMAGE_CHANNELS);

        assert(type == IL_FLOAT);
        assert(width == height);

        floatProbe.width = width;
        floatProbe.data = data;

        printf("%d  %d  %d\n", width, height, bytesPerPixel);
    }


    float coeffs[9][3] = {0} ;

    squareFloatImgBlob_reverseEndian(&floatProbe);


    prefilterAngleMap(&floatProbe, coeffs);

    printf("Filename            : %s\n", filename);
    printf("Filename (no ext)   : %s\n", getFilenameWithoutExt(filename));
    printf("Extension           : %s\n", getFilenameExt(filename));

    char* binaryFilename = (char*) malloc(strlen(filename) + 1 + strlen("shc") + 1);
    strcpy(binaryFilename, getFilenameWithoutExt(filename));
    strcat(binaryFilename, ".shc");

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


    printf("\n         Lighting Coefficients\n\n") ;
    printf("(l,m)       RED        GREEN     BLUE\n") ;

    printf("L_{0,0}   %9.6f %9.6f %9.6f\n",
       coeffs[0][0],coeffs[0][1],coeffs[0][2]) ;
    printf("L_{1,-1}  %9.6f %9.6f %9.6f\n",
       coeffs[1][0],coeffs[1][1],coeffs[1][2]) ;
    printf("L_{1,0}   %9.6f %9.6f %9.6f\n",
       coeffs[2][0],coeffs[2][1],coeffs[2][2]) ;
    printf("L_{1,1}   %9.6f %9.6f %9.6f\n",
       coeffs[3][0],coeffs[3][1],coeffs[3][2]) ;
    printf("L_{2,-2}  %9.6f %9.6f %9.6f\n",
       coeffs[4][0],coeffs[4][1],coeffs[4][2]) ;
    printf("L_{2,-1}  %9.6f %9.6f %9.6f\n",
       coeffs[5][0],coeffs[5][1],coeffs[5][2]) ;
    printf("L_{2,0}   %9.6f %9.6f %9.6f\n",
       coeffs[6][0],coeffs[6][1],coeffs[6][2]) ;
    printf("L_{2,1}   %9.6f %9.6f %9.6f\n",
       coeffs[7][0],coeffs[7][1],coeffs[7][2]) ;
    printf("L_{2,2}   %9.6f %9.6f %9.6f\n",
       coeffs[8][0],coeffs[8][1],coeffs[8][2]) ;

//    tomatrix(coeffs);

//    printMatricesToGlslDeclaration();


    return 0;
}


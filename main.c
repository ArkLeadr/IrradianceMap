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

struct squareByteImgBlob {
    int width;
    BYTE* data;
};
typedef struct squareByteImgBlob squareByteImgBlob;


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

void squareByteImgBlob_init(squareByteImgBlob* blob, int width) {
    blob->width = width;
    blob->data = malloc(width*width*3);
}

void squareFloatImgBlob_init(squareFloatImgBlob* blob, int width) {
    blob->width = width;
    blob->data = malloc(width*width*3*sizeof(float));
}

void convertByteImgToFloatImg(const squareByteImgBlob* inByteBlob, squareFloatImgBlob* outFloatBlob) {
    int width = inByteBlob->width;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int channel = 0; channel < 3; ++channel) {
                outFloatBlob->data[i * width * 3 + j * 3 + channel] = (float)inByteBlob->data[i * width * 3 + j * 3 + channel] / 255.f;
            }
        }
    }

}

void convertFloatImgToByteImgTonemap(const squareFloatImgBlob* inFloatBlob, squareByteImgBlob* outByteBlob) {
    int width = inFloatBlob->width;
    int pixelsCount = width * width;
    int totalSize = pixelsCount * 3; // 3 floats per pixel


    outByteBlob->width = inFloatBlob->width;
    outByteBlob->data = (BYTE*) malloc(totalSize);

    const double delta = 0.0001;

    float maxVal = 0.f;

    float lumMax[3] = {0.f};

    float avLum[3] = {0.f};

    float keyValue = 0.18f;

    float u, v, radius;

    int n = 0;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {

            v = (width/2.f - i)/(width/2.f);
            u = (j - width/2.f)/(width/2.f);
            radius = sqrtf(u*u + v*v) ;

            if (radius <= 1.0) {
                ++n;

                for (int channel = 0; channel < 3; ++channel) {
                    float currentLum = inFloatBlob->data[i * width * 3 + j * 3 + channel];

                    if (currentLum > lumMax[channel]) {
                        lumMax[channel] = currentLum;
                    }

                    avLum[channel] += log(currentLum + delta);
                }
            }

        }
    }

    avLum[0] = exp(avLum[0] / (float)n);
    avLum[1] = exp(avLum[1] / (float)n);
    avLum[2] = exp(avLum[2] / (float)n);



    float endMinLum[3] = {100000000.f, 100000000.f, 100000000.f};
    float endMaxLum[3] = {0.f};


//    n = 0.f;

//    for (int i = 0; i < width; ++i) {
//        for (int j = 0; j < width; ++j) {

//            v = (width/2.f - i)/(width/2.f);
//            u = (j - width/2.f)/(width/2.f);
//            radius = sqrtf(u*u + v*v) ;

//            if (radius <= 1.0) {
//                ++n;

//                for (int channel = 0; channel < 3; ++channel) {
//                    float floatComponent = inFloatBlob->data[i * width * 3 + j * 3 + channel];

//                    float currentLum = keyValue * floatComponent / avLum[channel];

//                    float lumWhite2 = lumMax[channel];

////                    currentLum = currentLum * (1.f + currentLum / lumWhite2) / (1.f + currentLum);

//                    currentLum = currentLum / (1.f + currentLum);

//                    if (currentLum > endMaxLum[channel]) endMaxLum[channel] = currentLum;
//                    if (currentLum < endMinLum[channel]) endMinLum[channel] = currentLum;
//                }
//            }

//        }
//    }


    for (int i = 0; i < totalSize; ++i) {
        int channel = i % 3;

        float floatComponent = inFloatBlob->data[i];

        float currentLum = keyValue * inFloatBlob->data[i] / avLum[channel];

//        currentLum = pow(currentLum, 2.2);

//        float lumWhite2 = lumMax[channel];

//        currentLum = currentLum * (1.f + currentLum / lumWhite2) / (1.f + currentLum);

        currentLum = currentLum / (1.f + currentLum);

//        currentLum = (currentLum - endMinLum[channel]) / (endMaxLum[channel] - endMinLum[channel]);

        outByteBlob->data[i] = (uint8_t) (currentLum * 255.f);

        if (channel == FI_RGBA_RED) {
            BYTE tmp = outByteBlob->data[i];
            outByteBlob->data[i] = outByteBlob->data[i - 2];
            outByteBlob->data[i - 2] = tmp;

//            outByteBlob->data[i] = min((int) outByteBlob->data[i] + 30, 255);
        }
        else {
//            outByteBlob->data[i] = 0;
        }
    }
}

void convertFloatImgToByteImg(const squareFloatImgBlob* inFloatBlob, squareByteImgBlob* outByteBlob) {
    int width = inFloatBlob->width;
    int pixelsCount = width * width;
    int totalSize = pixelsCount * 3; // 3 floats per pixel


    outByteBlob->width = inFloatBlob->width;
    outByteBlob->data = (BYTE*) malloc(totalSize);

    for (int i = 0; i < totalSize; ++i) {
        int channel = i % 3;

        float floatComponent = inFloatBlob->data[i];

        outByteBlob->data[i] = (uint8_t) (floatComponent);

        if (channel == FI_RGBA_RED) {
            BYTE tmp = outByteBlob->data[i];
            outByteBlob->data[i] = outByteBlob->data[i - 2];
            outByteBlob->data[i - 2] = tmp;
        }

    }
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

void swapByte(BYTE* lhs, BYTE* rhs) {
    BYTE tmp = *lhs;

    *lhs = *rhs;
    *rhs = tmp;
}

void squareByteImgBlob_loadFromDib(squareByteImgBlob* blob, FIBITMAP* dib) {
    squareByteImgBlob_init(blob, FreeImage_GetWidth(dib));

    FreeImage_ConvertToRawBits(blob->data, dib, FreeImage_GetWidth(dib) * 3, 24, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);

    // Swapping blue and red channel since FIF is BGR and masks seems not to work
    for (int i = 0; i < blob->width; ++i) {
        for (int j = 0; j < blob->width; ++j) {
            swapByte(&(blob->data[i * blob->width * 3 + j * 3 + 0]), &(blob->data[i * blob->width * 3 + j * 3 + 2]));
        }
    }
}

#define assignIfInf(var, minVar) if (var < minVar) minVar = var;
#define assignIfSup(var, maxVar) if (var > maxVar) maxVar = var;

#define between(inf, val, sup) (inf <= val && val <= sup)

#define approx(val, approximant) (between(approximant - 0.005, val, approximant + 0.005))

enum {
    POSX = 0,
    NEGX,
    POXY,
    NEGY,
    POSZ,
    NEGZ
};

const char* cubemap_filenames[6] = {
    "posx.jpg",
    "negx.jpg",
    "posy.jpg",
    "negy.jpg",
    "posz.jpg",
    "negz.jpg"
};

void stealOneQuarter(const squareFloatImgBlob* inFloatBlob, squareFloatImgBlob* quarter, int which) {
    double thetaStart;
    double phiStart;

    int width = inFloatBlob->width;

    switch(which) {
    case POSX:
        thetaStart = (-M_PI_4) + M_PI_2;
        phiStart = (-M_PI_4);
        break;
    case NEGX:
        thetaStart = (-M_PI_4) - M_PI_2;
        phiStart = (-M_PI_4);
        break;
    case POXY:
        thetaStart = (-M_PI_4);
        phiStart = (M_PI_4);
        break;
    case NEGY:
        thetaStart = (-M_PI_4);
        phiStart = (-M_PI_4 - M_PI_2);
        break;
    case POSZ:
        thetaStart = (-M_PI_4);
        phiStart = (-M_PI_4);
        break;
    case NEGZ:
        thetaStart = (-M_PI_4) + M_PI;
        phiStart = (-M_PI_4);
        break;
    }

    for (int i = 0 ; i < quarter->width ; i++) {
        for (int j = 0 ; j < quarter->width ; j++) {

            double r,theta,phi ;
            int fromI, fromJ;

            double x, y, z;

            theta = thetaStart + (M_PI_2 * ((double)j / (double)quarter->width));
            phi = phiStart + (M_PI_2 * ((double)i / (double)quarter->width));
            r = 1.0;

            z = r * cos(phi) * cos(theta);
            x = r * cos(phi) * sin(theta);
            y = r * sin(phi);

//            float multiplier = 1.f / x;

//            x *= multiplier;
//            y *= multiplier;
//            z *= multiplier;

            double norm = 1.0 / sqrt(x*x + y*y + z*z);

            double DDx = x * norm;
            double DDy = y * norm;
            double DDz = z * norm;

            double rr = 0.159154943 * acos(DDz) / sqrt(DDx*DDx + DDy*DDy);

            double sb_u = 0.5 + DDx * rr;
            double sb_v = 0.5 + DDy * rr;

            fromI = sb_v * width;
            fromJ = sb_u * width;

            for (int channel = 0; channel < 3; ++channel) {
                quarter->data[i * quarter->width * 3 + j * 3 + channel] = inFloatBlob->data[fromI * width * 3 + fromJ * 3 + channel];
            }
        }
    }

//    for (int i = 0 ; i < quarter->width ; i++) {
//        for (int j = 0 ; j < quarter->width ; j++) {

//            float u,v,r,theta,phi ;
//            int fromI, fromJ;

//            float x, y, z;

//            theta = (-M_PI_4) + (M_PI_2 * ((float)j / (float)quarter->width));
//            phi = (-M_PI_4) + (M_PI_2 * ((float)i / (float)quarter->width));
//            r = 1.f;

//            z = r * cosf(phi) * cosf(theta);
//            x = r * cosf(phi) * sinf(theta);
//            y = r * sinf(phi);

////            float multiplier = 1.f / x;

////            x *= multiplier;
////            y *= multiplier;
////            z *= multiplier;

//            float norm = 1.f / sqrtf(x*x + y*y + z*z);

//            float DDx = x * norm;
//            float DDy = y * norm;
//            float DDz = z * norm;

//            float rr = 0.159154943 * acosf(DDz) / sqrtf(DDx*DDx + DDy*DDy);

//            float sb_u = 0.5 + DDx * rr;
//            float sb_v = 0.5 + DDy * rr;

//            fromI = sb_v * width;
//            fromJ = sb_u * width;

//            for (int channel = 0; channel < 3; ++channel) {
//                quarter->data[i * quarter->width * 3 + j * 3 + channel] = inFloatBlob->data[fromI * width * 3 + fromJ * 3 + channel];
//            }
//        }
//    }

}


int main(void)
{
//    FIBITMAP* dib = FreeImage_Load(FIF_BMP, "grace_probe.bmp", 0);

//    squareByteImgBlob byteBmpBlob;
//    squareByteImgBlob_loadFromDib(&byteBmpBlob, dib);

//    squareFloatImgBlob floatBmpBlob;
//    squareFloatImgBlob_init(&floatBmpBlob, FreeImage_GetWidth(dib));

//    convertByteImgToFloatImg(&byteBmpBlob, &floatBmpBlob);

    squareFloatImgBlob floatProbe;

    loadFromFloatFormat("grace_probe.float", &floatProbe);

    float coeffs[9][3] = {0} ;

//    prefilterAngleMap(&floatProbe, coeffs);

//    prefilterAngleMap(&floatBmpBlob, coeffs);


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

//    printf("%p\n", bytes);

//    FIBITMAP* dibconv = FreeImage_ConvertFromRawBits(bytes, FreeImage_GetWidth(dib), FreeImage_GetWidth(dib), FreeImage_GetWidth(dib) * 3,
//                                                     24,
//                                                     FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);

//    FreeImage_Save(FIF_BMP, dibconv, "test.bmp", 0);

    squareByteImgBlob byteProbe;
    squareFloatImgBlob quarter;

    squareFloatImgBlob_init(&quarter, 1024);

    for (int i = 0; i < 6; ++i) {
        stealOneQuarter(&floatProbe, &quarter, i);

        convertFloatImgToByteImgTonemap(&quarter, &byteProbe);

        FIBITMAP* dibsave = FreeImage_ConvertFromRawBits(byteProbe.data, byteProbe.width, byteProbe.width, byteProbe.width * 3, 24,
                                                    FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);

        FreeImage_Save(FIF_JPEG, dibsave, cubemap_filenames[i], 0);
    }

    return 0;
}


//To compile and run:
// On os x:
// g++ cleanLinearCombinationInterpolation.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o cleanLinearCombinationInterpolation -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./cleanLinearCombinationInterpolation
//
// On linux:
// g++ cleanLinearCombinationInterpolation.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o cleanLinearCombinationInterpolation -std=c++11
//
// LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./cleanLinearCombinationInterpolation

//uses tuples (linearCombinationKernel does not)
//single file (not modularized like linearCombinationKernel1)
//this kernel is a spatially varying linear combination of guassians 
//that uses tuples for fast evaluation

#define STANDALONE
#ifndef STANDALONE
#include "lsst/afw/image.h"
namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
#endif

#include <stdio.h>
#include "Halide.h"
#include <bitset>
#include "clock.h"
using namespace std;
using namespace Halide;

using Halide::Image;

int main(int argc, char *argv[]) {

#ifndef STANDALONE
    auto im = afwImage::MaskedImage<float>("./images/calexp-004207-g3-0123.fits");
    int width = im.getWidth(), height = im.getHeight();
#else
    int width = 2048, height = 1489;
    printf("[no load]");
#endif
    printf("Loaded: %d x %d\n", width, height);

    //store image data in img_var(x, y, 0) and variance data in img_var(x, y, 1)
    Image<float> image(width, height);
    Image<float> variance(width, height);
    Image<uint16_t> mask(width, height);

#ifndef STANDALONE
    //Read image in
    for (int y = 0; y < height; y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = im.x_at(0, y);
        for (int x = 0; x < width; x++){
            image(x, y) = (*inPtr).image();
            variance(x, y) = (*inPtr).variance();
            mask(x, y) = (*inPtr).mask();
            inPtr++;
        }
    }
#endif

    //Kernel has dimensions (boundingBox*2 + 1) x (boundingBox*2 + 1)
    int boundingBox = 2; 
    //Interpolate the kernel over an interpDistxinterpDist grid where the
    //exact kernel is computed
    int interpDist = 10;
    Var x, y, i, j, y0, yi;

    //Five 3rd degree polynomials which will be used as spatially varying
    //coefficients in the linear combination of the five gaussian basis kernels
    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.002f*x + 0.003f*y + 0.4f*x*x + 0.5f*x*y
                     + 0.6f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
                     + 0.00011f*y*y*y;

    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 1.1f + 1.002f*x + 1.003f*y + 1.4f*x*x + 1.5f*x*y
                     + 1.6f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
                     + 1.00011f*y*y*y;

    Func polynomial3 ("polynomial3");
    polynomial3(x, y) = 2.1f + 2.002f*x + 2.003f*y + 2.4f*x*x + 2.5f*x*y
                     + 2.6f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
                     + 2.00011f*y*y*y;

    Func polynomial4 ("polynomial4");
    polynomial4(x, y) = 3.1f + 3.002f*x + 3.003f*y + 3.4f*x*x + 3.5f*x*y
                     + 3.6f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
                     + 3.00011f*y*y*y;

    Func polynomial5 ("polynomial5");
    polynomial5(x, y) = 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
                     + 4.00011f*y*y*y;

    //5 Guassian basis kernels
    float pi = 3.14159265359f;
    //Kernel #1
    Func kernel1;
    float sigmaX1 = 2.0f;
    float sigmaY1 = 2.0f;
    float theta1 = 0.0f; //rotation of sigmaX axis
    kernel1(x, y) = exp(-((x*cos(theta1) + y*sin(theta1))*(x*cos(theta1) + y*sin(theta1)))
                    /(2*sigmaX1*sigmaX1)
                    -((y*cos(theta1) - x*sin(theta1))*(y*cos(theta1) - x*sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (2.0f*pi*sigmaX1*sigmaY1);

    //Kernel #2
    Func kernel2;
    float sigmaX2 = 0.5f;
    float sigmaY2 = 4.0f;
    float theta2 = 0.0f; //rotation of sigmaX axis
    kernel2(x, y) = exp(-((x*cos(theta2) + y*sin(theta2))*(x*cos(theta2) + y*sin(theta2)))
                    /(2*sigmaX2*sigmaX2)
                    -((y*cos(theta2) - x*sin(theta2))*(y*cos(theta2) - x*sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (2.0f*pi*sigmaX2*sigmaY2);

    //Kernel #3
    Func kernel3;
    float sigmaX3 = 0.5f;
    float sigmaY3 = 4.0f;
    float theta3 = 3.14159f/4; //rotation of sigmaX axis
    kernel3(x, y) = exp(-((x*cos(theta3) + y*sin(theta3))*(x*cos(theta3) + y*sin(theta3)))
                    /(2*sigmaX3*sigmaX3)
                    -((y*cos(theta3) - x*sin(theta3))*(y*cos(theta3) - x*sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (2.0f*pi*sigmaX3*sigmaY3);
    //Kernel #4
    Func kernel4;
    float sigmaX4 = 0.5f;
    float sigmaY4 = 4.0f;
    float theta4 = 3.14159f/2; //rotation of sigmaX axis
    kernel4(x, y) = exp(-((x*cos(theta4) + y*sin(theta4))*(x*cos(theta4) + y*sin(theta4)))
                    /(2*sigmaX4*sigmaX4)
                    -((y*cos(theta4) - x*sin(theta4))*(y*cos(theta4) - x*sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (2.0f*pi*sigmaX4*sigmaY4);


    //Kernel #5
    Func kernel5;
    float sigmaX5 = 4.0f;
    float sigmaY5 = 4.0f;
    float theta5 = 0.0; //rotation of sigmaX axis
    kernel5(x, y) = (exp(-((x*cos(theta5) +y*sin(theta5))*(x*cos(theta5) +y*sin(theta5)))
                    /(2*sigmaX5*sigmaX5)))
                    *(exp(-((y*cos(theta5) - x*sin(theta5))*(y*cos(theta5) - x*sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (2.0f*pi*sigmaX5*sigmaY5));

    //compute the normalized linear combination of basis kernels 
    //(so the kernel sums to 1)
    Func normalizedTotalKernel;
    Expr norm = 0.0f;

    for(int ii = -boundingBox; ii <= boundingBox; ii++){
        for(int jj = -boundingBox; jj <= boundingBox; jj++){
            norm += (polynomial1(x, y)*kernel1(ii, jj) +
                polynomial2(x, y)*kernel2(ii, jj) + polynomial3(x, y)*kernel3(ii, jj) + 
                polynomial4(x, y)*kernel4(ii, jj) + polynomial5(x, y)*kernel5(ii, jj));
        }
    }

    normalizedTotalKernel(x, y, i, j) = (polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j))/norm;

    //This is a compressed version of the original kernel that can be evaluated
    //compute_root() without computing the kernel at points that will not
    //actually be used for interpolation.
    Func compressedKernel;
    compressedKernel(x, y, i, j) = normalizedTotalKernel(x*interpDist, y*interpDist, i, j);

    //Create a kernel that is the 2d linear interpolation of the compressed kernel
    Func interpKernel;
    Expr xS = x%interpDist;
    Expr xG = interpDist - xS;
    Expr yS = y%interpDist;
    Expr yG = interpDist - yS;

    Expr xInterpS = ((xG * compressedKernel((x - xS)/interpDist, (y - yS)/interpDist, i, j)) +
                     (xS * compressedKernel((x + xG)/interpDist, (y - yS)/interpDist, i, j))) / interpDist;

    Expr xInterpG = ((xG * compressedKernel((x - xS)/interpDist, (y + yG)/interpDist, i, j)) +
                     (xS * compressedKernel((x + xG)/interpDist, (y + yG)/interpDist, i, j))) / interpDist;

    interpKernel(x, y, i, j) = ((yG * xInterpS) + (yS * xInterpG))/interpDist;





    //Compute output image plane
    Func image_bounded ("image_bounded");
    image_bounded = BoundaryConditions::repeat_edge(image);
    Expr imageOutHelp = 0.0f;

    //Compute output variance plane
    Func variance_bounded ("variance_bounded");
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    Expr varianceOutHelp = 0.0f;

    //Compute output mask plane
    Func mask_bounded ("mask_bounded");
    mask_bounded = BoundaryConditions::repeat_edge(mask);
    Expr maskOutHelp = cast<uint16_t>(0);

    Expr curKernelVal;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            imageOutHelp += image_bounded(x + i, y + j)*interpKernel(x, y, i, j);
            varianceOutHelp += variance_bounded(x + i, y + j)*interpKernel(x, y, i, j)
                                    *interpKernel(x, y, i, j);

            maskOutHelp = select(interpKernel(x, y, i, j) == 0.0f, maskOutHelp,
                                    maskOutHelp | mask_bounded(x + i, y + j));
        }
    }


    //Evaluate image, mask, and variance planes concurrently using a tuple
    Func combined_output ("combined_output");
    combined_output(x, y) = Tuple(imageOutHelp, varianceOutHelp, maskOutHelp);


    //Schedule:
    // Split the y coordinate of the output into strips of 32 scanlines:
    combined_output.split(y, y0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y0);
    // Vectorize across x by a factor of eight.
    combined_output.vectorize(x, 8);

    compressedKernel.compute_root();

    //Create output images for the image, variance, and mask planes
    Image<float> image_output(image.width(), image.height());
    Image<float> variance_output(variance.width(), variance.height());
    Image<uint16_t> mask_output(mask.width(), mask.height());

    //Compute all three planes simultaneously using a tuple
    //Evaluate once before benchmarking to force Halide to compile
    Realization r = combined_output.realize(image.width(), image.height());
    //Pull the three output planes out of the tuple
    image_output = r[0];
    variance_output = r[1];
    mask_output = r[2];

    // Benchmark the pipeline.
    double mean = 0;
    double min;
    double max;
    int numberOfRuns = 5;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        r = combined_output.realize(image.width(), image.height());
        double t2 = current_time();
        double curTime = (t2-t1);
        mean += curTime;
        if(i == 0){
            min = curTime;
            max = curTime;
        }
        else{
            if(curTime < min)
                min = curTime;
            if(curTime > max)
                max = curTime;
        }
    }
    mean = mean/numberOfRuns;

    std::cout << "Mean Time: " << mean << ", Min = " <<
    min << ", Max = " << max << ", with " << numberOfRuns <<
    " runs" << '\n';

#ifndef STANDALONE
    //write image out
    auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
    for (int y = 0; y < imOut.getHeight(); y++) {
    	afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);

        for (int x = 0; x < imOut.getWidth(); x++){
        	afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
            curPixel(image_output(x, y), mask_output(x, y), variance_output(x, y));
        	(*inPtr) = curPixel;
        	inPtr++;

        }
    }

	imOut.writeFits("./halideCleanLinearCombinationInterpolation5x5.fits");
#endif

}


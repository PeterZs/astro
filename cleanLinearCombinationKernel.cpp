//To compile and run:
// On os x:
// g++ cleanLinearCombinationKernel.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o cleanLinearCombinationKernel -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./cleanLinearCombinationKernel
//
// On linux:
// g++ cleanLinearCombinationKernel.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o cleanLinearCombinationKernel -std=c++11
//
// LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./cleanLinearCombinationKernel

//uses tuples (linearCombinationKernel does not)
//single file (not modularized like linearCombinationKernel1)
//this kernel is a spatially varying linear combination of guassians 
//that uses tuples for fast evaluation

//#define STANDALONE
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

#ifdef NAN
    cout << NAN << endl;
    float test = NAN;
    cout << test << endl;
#else
    cout << "nan not defined" << endl;
#endif


#ifndef STANDALONE
    auto im = afwImage::MaskedImage<float>("./images/calexp-004207-g3-0123.fits");
    int width = im.getWidth(), height = im.getHeight();
#else
    int width = 2048, height = 1489;
    printf("[no load]");
#endif
    printf("Loaded: %d x %d\n", width, height);

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


    //Compute output image plane
    Func image_bounded ("image_bounded");
    image_bounded = BoundaryConditions::repeat_edge(image);
    Expr blur_image_help = 0.0f;
    Expr norm = 0.0f;

    //Compute output variance plane
    Func variance_bounded ("variance_bounded");
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    Func blurVariance ("blurVariance");
    Expr blur_variance_help = 0.0f;

    //Compute output mask plane
    Func mask_bounded ("mask_bounded");
    mask_bounded = BoundaryConditions::repeat_edge(mask);
    Func maskOut ("maskOut");
    Expr maskOutHelp = cast<uint16_t>(0);

//Slow using reductions below:
//*********************************************************************************

    Expr curKernelVal = (polynomial1(x, y)*kernel1(i, j) +
        polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
        polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j));

    Func cur_Kernel_Val;
    cur_Kernel_Val(x, y, i, j) = curKernelVal;

    blur_image_help = image_bounded(x + i, y + j)*curKernelVal; 
    blur_variance_help = variance_bounded(x + i, y + j)*curKernelVal*curKernelVal; 
    maskOutHelp = select(curKernelVal == 0.0f, maskOutHelp,
                        maskOutHelp | mask_bounded(x + i, y + j));


    Func blur_image_help_func;
    blur_image_help_func(x, y, i, j) = blur_image_help;
    Func blur_variance_help_func;
    blur_variance_help_func(x, y, i, j) = blur_variance_help;
    Func blur_mask_help_func;
    blur_mask_help_func(x, y, i, j) = maskOutHelp;

    RDom r(-2, 5, -2, 5);
    //Evaluate image, mask, and variance planes concurrently using a tuple
    Func combined_output ("combined_output");

    Func norm_func;
//    norm_func(x, y) = 0.0f;
//    norm_func(x, y) += cur_Kernel_Val(x, y, r.x, r.y);
    norm_func(x, y) = sum(cur_Kernel_Val(x, y, r.x, r.y));

    //set the image edges
    //image edge should be NAN, but this produces errors 
    Func setEdge;
    setEdge(x, y) = x < boundingBox || y < boundingBox ||
                     x > (width - 1 - boundingBox) || y > (height - 1 - boundingBox);


    Func image_output_func;
//    image_output_func(x, y) = 0.0f;
//    image_output_func(x, y) += blur_image_help_func(x, y, r.x, r.y);
    image_output_func(x, y) = sum(blur_image_help_func(x, y, r.x, r.y));
    image_output_func(x, y) = image_output_func(x, y) / norm_func(x, y);
    image_output_func(x, y) = select(setEdge(x, y), INFINITY, image_output_func(x, y)); 

    Func var_output_func;
//    var_output_func(x, y) = 0.0f;
//    var_output_func(x, y) += blur_variance_help_func(x, y, r.x, r.y);
    var_output_func(x, y) = sum(blur_variance_help_func(x, y, r.x, r.y));
    var_output_func(x, y) = var_output_func(x, y) / (norm_func(x, y) * norm_func(x, y));
    var_output_func(x, y) = select(setEdge(x, y), INFINITY, var_output_func(x, y)); 

    Func mask_output_func;
//    mask_output_func(x, y) = cast<uint16_t>(0);
//    mask_output_func(x, y) += blur_mask_help_func(x, y, r.x, r.y);
    mask_output_func(x, y) = sum(blur_mask_help_func(x, y, r.x, r.y));
    mask_output_func(x, y) = select(setEdge(x, y), 16, mask_output_func(x, y)); 

    combined_output(x, y) = Tuple(image_output_func(x, y), var_output_func(x, y)
                                    , mask_output_func(x, y));

//*********************************************************************************











//Fast without reductions below:
//*********************************************************************************

//    Expr curKernelVal;
//    for(int i = -boundingBox; i <= boundingBox; i++){
//        for(int j = -boundingBox; j <= boundingBox; j++){
//
//            curKernelVal = (polynomial1(x, y)*kernel1(i, j) +
//                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
//                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j));
//            blur_image_help += image_bounded(x + i, y + j)*curKernelVal; 
//            blur_variance_help += variance_bounded(x + i, y + j)*curKernelVal*curKernelVal; 
//            maskOutHelp = select(curKernelVal == 0.0f, maskOutHelp,
//                                maskOutHelp | mask_bounded(x + i, y + j));
//
//
//            norm += curKernelVal;
//        }
//    }
//    blur_image_help = blur_image_help/norm;
//    blur_variance_help = blur_variance_help/(norm*norm);
//
//    //set the image edges
//    //image edge should be NAN, but this produces errors 
//    Expr setEdge = x < boundingBox || y < boundingBox ||
//                     x > (width - 1 - boundingBox) || y > (height - 1 - boundingBox);
//    blur_image_help = select(setEdge, INFINITY, blur_image_help); 
//    blur_variance_help = select(setEdge, INFINITY, blur_variance_help);
//    maskOutHelp = select(setEdge, 16, maskOutHelp);
//
//    //Evaluate image, mask, and variance planes concurrently using a tuple
//    Func combined_output ("combined_output");
//    combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, maskOutHelp);
//*********************************************************************************

    //set the image edges
//    for(int i = 0; i <= boundingBox; i++){
//        combined_output(i, y)[0] = NAN;
//        combined_output(width-1-i, y)[0] = NAN;
//        combined_output(x, i)[0] = NAN;
//        combined_output(x, height-1-i)[0] = NAN;
//
//        combined_output(i, y)[1] = INFINITY;
//        combined_output(width-1-i, y)[1] = INFINITY;
//        combined_output(x, i)[1] = INFINITY;
//        combined_output(x, height-1-i)[1] = INFINITY;
//
//        combined_output(i, y)[2] = 16.0f;
//        combined_output(width-1-i, y)[2] = 16.0f;
//        combined_output(x, i)[2] = 16.0f;
//        combined_output(x, height-1-i)[2] = 16.0f;
//
//    }


    // Split the y coordinate of the output into strips of 32 scanlines:
    combined_output.split(y, y0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y0);
    // Vectorize across x by a factor of eight.
    combined_output.vectorize(x, 8);

    //cur_Kernel_Val.compute_at(combined_output, x);
//    polynomial5.compute_root();
//    polynomial4.compute_root();
//    polynomial3.compute_root();
//    polynomial2.compute_root();
//    polynomial1.compute_root();

    //Create output images for the image, variance, and mask planes
    Image<float> image_output(image.width(), image.height());
    Image<float> variance_output(variance.width(), variance.height());
    Image<uint16_t> mask_output(mask.width(), mask.height());

    //Compute all three planes simultaneously using a tuple
    //Evaluate once before benchmarking to force Halide to compile
    Realization rOut = combined_output.realize(image.width(), image.height());
    //Pull the three output planes out of the tuple
    image_output = rOut[0];
    variance_output = rOut[1];
    mask_output = rOut[2];

    // Benchmark the pipeline.
    double mean = 0;
    double min;
    double max;
    int numberOfRuns = 5;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        rOut = combined_output.realize(image.width(), image.height());
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
//    //write image out
//    auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
//    for (int y = 0; y < imOut.getHeight(); y++) {
//    	afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);
//
//        for (int x = 0; x < imOut.getWidth(); x++){
//        	afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
//            curPixel(image_output(x, y), mask_output(x, y), variance_output(x, y));
//        	(*inPtr) = curPixel;
//        	inPtr++;
//
//        }
//    }
//
//	imOut.writeFits("./halideCleanLinearCombination5x5.fits");
        //write image out
    auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
    for (int y = 0; y < imOut.getHeight(); y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);

        for (int x = 0; x < imOut.getWidth(); x++){
            afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
            curPixel(image_output(x, y), mask(x, y), variance(x, y));
            (*inPtr) = curPixel;
            inPtr++;

        }
    }

    auto varOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
    for (int y = 0; y < imOut.getHeight(); y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = varOut.x_at(0, y);

        for (int x = 0; x < imOut.getWidth(); x++){
            afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
            curPixel(variance_output(x, y), mask(x, y), variance(x, y));
            (*inPtr) = curPixel;
            inPtr++;

        }
    }

    auto maskOutPlane = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
    for (int y = 0; y < imOut.getHeight(); y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = maskOutPlane.x_at(0, y);

        for (int x = 0; x < imOut.getWidth(); x++){
            afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
            curPixel(mask_output(x, y), mask(x, y), variance(x, y));
            (*inPtr) = curPixel;
            inPtr++;

        }
    }

    imOut.writeFits("./halideLinComboImage5x5.fits");
    varOut.writeFits("./halideLinComboVar5x5.fits");
    maskOutPlane.writeFits("./halideLinComboMask5x5.fits");
#endif

}


// On os x:
// g++ analyticKernelInterpolate3.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o analyticKernelInterpolate3 -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./analyticKernelInterpolate3
//
// On linux:
// g++ analyticKernelInterpolate3.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o analyticKernelInterpolate3 -std=c++11
//
//LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./analyticKernelInterpolate3


//this code performs interpolation on the Kernel itself as opposed to the underlying polynomial
//Image, mask, and variance are stored in separate halide images, so that when splitting/tiling updates there are no boundary condition
//violations.  Still need to figure out how to combine kernel computation for three planes.

#include "lsst/afw/image.h"
#include <stdio.h>

#include "Halide.h"
#include <bitset>
#include "clock.h"
using namespace std;
using namespace Halide;

using Halide::Image;

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;


int main(int argc, char *argv[]) {

    //This image has 3 planes (image, mask, variance) and dimensions 2048 x 1489
    auto im = afwImage::MaskedImage<float>("./images/calexp-004207-g3-0123.fits");
    printf("Loaded: %d x %d\n", im.getWidth(), im.getHeight());

    //store image data in img_var(x, y, 0) and variance data in img_var(x, y, 1)
    Image<float> image(im.getWidth(), im.getHeight());
    Image<float> variance(im.getWidth(), im.getHeight());
    Image<uint16_t> mask(im.getWidth(), im.getHeight());

    //Read image in
    for (int y = 0; y < im.getHeight(); y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = im.x_at(0, y);
        for (int x = 0; x < im.getWidth(); x++){
            image(x, y) = (*inPtr).image();
            variance(x, y) = (*inPtr).variance();
            mask(x, y) = (*inPtr).mask();
            inPtr++;
        }
    }

    Func image_bounded ("image_bounded");
    image_bounded = BoundaryConditions::repeat_edge(image);
    Func variance_bounded ("variance_bounded");
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    Func mask_bounded ("mask_bounded");
    mask_bounded = BoundaryConditions::repeat_edge(mask);

    int boundingBox = 2; 
    int interpConst = 10; //interpolate polynomial over interpConst x interpConst box
    Var x, y, i, j, i_v, y0, yi, x_outer, y_outer, x_inner, y_inner, fused, tile_index;

    //compute output image and variance
    Func polynomial ("polynomial");
    polynomial(x, y) = 0.1f + 0.0f*x + 0.0019476158495634653f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.0f*x + 0.0019476158495634653f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 0.1f + 0.0f*x + 0.0019476158495634653f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;



    //The kernel is spatially variant for different positions (x, y) in the image
    //kernel(x, y, i, j) is the value of the kernel at position (x, y) in the image
    //and position (i, j) inside the kernel
    Func kernel ("kernel");
    kernel(x, y, i, j) = (exp(-((i*cos(polynomial(x, y)) +j*sin(polynomial(x, y)))
                    *(i*cos(polynomial(x, y)) +j*sin(polynomial(x, y))))
                    /(2.0f*polynomial(x, y)*polynomial(x, y))) / (sqrtf(2.0f*M_PI)*polynomial(x, y)))
                    *(exp(-((j*cos(polynomial(x, y)) - i*sin(polynomial(x, y)))
                    *(j*cos(polynomial(x, y)) - i*sin(polynomial(x, y))))
                    /(2.0f*polynomial(x, y)*polynomial(x, y))) / (sqrtf(2.0f*M_PI)*polynomial(x, y)));

    //for experimenting with optimizations
    Func kernel1 ("kernel1");
    kernel1(x, y, i, j) = (exp(-((i*cos(polynomial1(x, y)) +j*sin(polynomial1(x, y)))
                    *(i*cos(polynomial1(x, y)) +j*sin(polynomial1(x, y))))
                    /(2.0f*polynomial1(x, y)*polynomial1(x, y))) / (sqrtf(2.0f*M_PI)*polynomial1(x, y)))
                    *(exp(-((j*cos(polynomial1(x, y)) - i*sin(polynomial1(x, y)))
                    *(j*cos(polynomial1(x, y)) - i*sin(polynomial1(x, y))))
                    /(2.0f*polynomial1(x, y)*polynomial1(x, y))) / (sqrtf(2.0f*M_PI)*polynomial1(x, y)));
    
    Func kernel2 ("kernel2");
    kernel2(x, y, i, j) = (exp(-((i*cos(polynomial2(x, y)) +j*sin(polynomial2(x, y)))
                    *(i*cos(polynomial2(x, y)) +j*sin(polynomial2(x, y))))
                    /(2.0f*polynomial2(x, y)*polynomial2(x, y))) / (sqrtf(2.0f*M_PI)*polynomial2(x, y)))
                    *(exp(-((j*cos(polynomial2(x, y)) - i*sin(polynomial2(x, y)))
                    *(j*cos(polynomial2(x, y)) - i*sin(polynomial2(x, y))))
                    /(2.0f*polynomial2(x, y)*polynomial2(x, y))) / (sqrtf(2.0f*M_PI)*polynomial2(x, y)));
//******************************VALID ORIGINAL
/*    Func interpKernel;
    Expr xS = x%interpConst;
    Expr xG = interpConst - xS;
    Expr yS = y%interpConst;
    Expr yG = interpConst - yS;

    Expr xInterpS = ((xG * kernel1(x - xS, y - yS, i, j)) + (xS * kernel1(x + xG, y - yS, i, j)))/interpConst;
    Expr xInterpG = ((xG * kernel1(x - xS, y + yG, i, j)) + (xS * kernel1(x + xG, y + yG, i, j)))/interpConst;
*/
//******************************DONE VALID ORIGINAL


//**************************COMPRESSED GRID (incorrect kernel, but speed should be ok)
    Func interpKernel;
    Expr xS = x%interpConst;
    Expr xG = interpConst - xS;
    Expr yS = y%interpConst;
    Expr yG = interpConst - yS;

    Expr xInterpS = ((xG * kernel1((x - xS)/interpConst, (y - yS)/interpConst, i, j)) + (xS * kernel1((x + xG)/interpConst, (y - yS)/interpConst, i, j)))/interpConst;
    Expr xInterpG = ((xG * kernel1((x - xS)/interpConst, (y + yG)/interpConst, i, j)) + (xS * kernel1((x + xG)/interpConst, (y + yG)/interpConst, i, j)))/interpConst;

//**************************Done COMPRESSED GRID


    interpKernel(x, y, i, j) = ((yG * xInterpS) + (yS * xInterpG))/interpConst;

    //compute Image output
    Func blurImage ("blurImage");
    Expr blur_image_help = 0.0f;
    Expr norm1 = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            blur_image_help += image_bounded(x + i, y + j) * interpKernel(x, y, i, j); 
            norm1 += interpKernel(x, y, i, j);
//            blur_image_help += image_bounded(x + i, y + j) * kernel(x, y, i, j); 
//            norm1 += kernel(x, y, i, j);

        }
    }

    blur_image_help = blur_image_help/norm1;

    blurImage(x, y) = blur_image_help;


    //compute Variance output
    Func blurVariance ("blurVariance");
    Expr blur_variance_help = 0.0f;
    Expr norm2 = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            blur_variance_help += variance_bounded(x + i, y + j) * kernel(x, y, i, j) * kernel(x, y, i, j); 
            norm2 += kernel(x, y, i, j);
        }
    }
//    blur_variance_help = blur_variance_help/(norm(x,y)*norm(x,y));
    blur_variance_help = blur_variance_help/(norm2*norm2);
    blurVariance(x, y) = blur_variance_help;

    //Compute mask
    Func maskOut ("maskOut");
    Expr maskOutHelp = 0;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            maskOutHelp = select(kernel(5, 5, i, j) == 0.0f, maskOutHelp, maskOutHelp | mask_bounded(x + i, y + j));
//            maskOutHelp = maskOutHelp | mask_bounded(x + i, y + j);      
        }
    }

    maskOut(x, y) = maskOutHelp;

    //Schedule

//        kernel1.reorder(i, j, x, y);
//        interpKernel.reorder(i, j, x, y);
//        kernel.reorder(i, j, x, y);
//        polynomial1.compute_root();
 //       kernel1.compute_root();
//        kernel1.compute_root();
//        kernel1.compute_at(blurImage, x).store_at(blurImage, y0);
//        kernel1.compute_at(blurImage, x_inner).store_at(blurImage, tile_index);

    
    
        blurImage.split(y, y0, yi, 10);
        blurImage.parallel(y0);
        blurImage.vectorize(x, 8);

/*        blurImage
        .tile(x, y, x_outer, y_outer, x_inner, y_inner, 20, 20)
            .fuse(x_outer, y_outer, tile_index)
            .parallel(tile_index)
//            .vectorize(x_inner, 4);
            .fuse(x_inner, y_inner, fused)
            .vectorize(fused, 4)
            .unroll(fused, 25);
*/

        blurVariance.split(y, y0, yi, 10);
        blurVariance.parallel(y0);
        blurVariance.vectorize(x, 8);


    maskOut.split(y, y0, yi, 4);
    maskOut.parallel(y0);
    maskOut.vectorize(x, 8);

//    kernel1.trace_stores();
//    blurImage.trace_stores();
    blurImage.print_loop_nest();
    blurVariance.print_loop_nest();
    maskOut.print_loop_nest();

    kernel1.trace_realizations();


    // Print out pseudocode for the pipeline.
//    blurImage.compile_to_lowered_stmt("./InterpolateSchedules/imageSchedule.html", {image}, HTML);
//    blurVariance.compile_to_lowered_stmt("analyticKernelInterpolate3BlurVariance.html", {variance}, HTML);
//    maskOut.compile_to_lowered_stmt("analyticKernelInterpolate3BlurMask.html", {mask}, HTML);
//    blurImage.compile_to_c("analyticKernelInterpolate3_C_Code.cpp", std::vector<Argument>(), "analyticKernelInterpolate3_C_Code");
//    blurVariance.compile_to_lowered_stmt("blur.html", {variance}, HTML);


    // Benchmark the pipeline.
    Image<float> image_output(image.width(), image.height());
    blurImage.realize(image_output);

    Image<float> variance_output(variance.width(), variance.height());
    blurVariance.realize(variance_output);

    Image<int32_t> mask_output(mask.width(), mask.height());
    maskOut.realize(mask_output);

    double average = 0;
    double min;
    double max;
    double imgTime;
    double varTime;
    double maskTime;
    int numberOfRuns = 1;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        blurImage.realize(image_output);
        double t2 = current_time();
        blurVariance.realize(variance_output);
        double t3 = current_time();
        maskOut.realize(mask_output);
        double t4 = current_time();
        double curTime = (t4-t1);
        average += curTime;
        if(i == 0){
            min = curTime;
            max = curTime;
            imgTime = t2-t1;
            varTime = t3-t2;
            maskTime = t4-t3;
        }
        else{
            if(curTime < min){
                min = curTime;
                imgTime = t2-t1;
                varTime = t3-t2;
                maskTime = t4-t3;
            }
            if(curTime > max)
                max = curTime;
        }
    }
    average = average/numberOfRuns;
    std::cout << "Average Time: " << average << ", Min = " <<
    min << ", Max = " << max << ", with " << numberOfRuns <<
    " runs" << '\n';
    cout << "For fastest run total time = " << min << ", imgTime = " << imgTime << ", varTime = " << varTime << 
    " maskTime = " << maskTime << endl;

//    blur_mask.realize(mask_output);

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

    imOut.writeFits("./images/halideAnalyticKernelInterpolate3.fits");
}



//    polynomial2(x, y) = 0.05f + 0.049f; slow
//    polynomial2(x, y) = 0.05f + 0.046f; not slow
//    polynomial2(x, y) = 0.05f + 0.04f; not slow

//    polynomial2(x, y) = 0.05f + 0.055f; slow
//    polynomial2(x, y) = 0.05f + 0.056f;
//    polynomial2(x, y) = 0.05f + 0.0565f; varies between slow/not slow up and down
//    polynomial2(x, y) = 0.05f + 0.057f;
//    polynomial2(x, y) = 0.05f + 0.058f; not slow



//*****************FAKE TESTING
/*    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
//            blur_image_help += image_bounded(x + i, y + j) * kernel2(x%10 + 1, y%10 + 1, i, j); 
//            norm1 += kernel2(x%10 + 1, y%10 + 1, i, j);

            blur_image_help += image_bounded(x + i, y + j) * kernel(348, 348, i, j); 
            norm1 += kernel(348, 348, i, j);

//            blur_image_help += image_bounded(x + i, y + j) * kernel(x, y, i, j); 
//            norm1 += kernel(x, y, i, j);

        }
    }
*/
//*****************DONE FAKE TESTING


 //  blur.reorder(i_v, x, y);


//    kernel1.compute_at(blurImage, x);
//    kernel1.vectorize(x, 8);
//    kernel1.split(y, y0, yi, 4);
//    kernel1.parallel(y0);


    //best schedule found:

/*
    // Split the y coordinate of the consumer into strips:
//    blurImage.split(y, y0, yi, 4);
    blurVariance.split(y, y0, yi, 4);
    // Compute the strips using a thread pool and a task queue.
//    blurImage.parallel(y0);
    blurVariance.parallel(y0);
    // Vectorize across x.
//    blurImage.vectorize(x, 8);
    blurVariance.vectorize(x, 8);

    kernel1.reorder(i, j, x, y);
    interpKernel.reorder(i, j, x, y);
    kernel.reorder(i, j, x, y);


//    kernel2.reorder(i, j, x, y);
//    kernel2.compute_root();

//    kernel1.compute_root();
//    kernel1.vectorize(x, 8);
//    kernel1.split(y, y0, yi, 4);
//    kernel1.parallel(y0);
//    kernel1.compute_at(blurImage, x).store_at(blurImage, y0);
//    interpKernel.compute_at(blurImage, x);



//    blurImage.split(x, x_outer, x_inner, 10);
//    blurImage.unroll(x_inner);



    blurImage.tile(x, y, x_outer, y_outer, x_inner, y_inner, 100, 10)
        .fuse(x_outer, y_outer, tile_index)
        .parallel(tile_index)
        .fuse(x_inner, y_inner, fused)
        .unroll(fused, 10);
//        .vectorize(x_inner, 8);
//        .fuse(x_inner, y_inner, fused)
//        .unroll(fused, 10)
//        .parallel(y_outer);
//        .vectorize(x_inner, 8);

//    kernel1.compute_at(blurImage, fused)
//        .store_at(blurImage, fused)
//        .vectorize(x, 4);
//    kernel1.compute_at(blurImage, y_inner)
//        .store_at(blurImage, y_inner)
//        .vectorize(x, 8);

//
//    kernel1.compute_at(blurImage, y_inner)
//        .store_at(blurImage, y_inner)
//        .vectorize(x, 4);


//    polynomial1.compute_at(blurImage, x).vectorize(x, 8);
//    kernel1.compute_at(blurImage, x).vectorize(x, 8);
*/



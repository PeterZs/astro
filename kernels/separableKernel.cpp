//for testing without output image 
//compile using: $make blurFast.standalone
//run using: ./blurFast.standalone

//for testing with output image
// On os x:
// g++ blurFast.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o blurFast -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./blurFast
//
// On linux:
// g++ blurFast.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o blurFast -std=c++11
//
//LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./blurFast



//jit compilation is hanging using vectorization and a straightforward loop 
//expression of blur_y (or blur_x).  If the kernel terms in the expression of 
//blur_y are reordered or vectorization is removed, the code compiles.


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
    auto im = afwImage::MaskedImage<float>("../calexp-004207-g3-0123.fits");
    int width = im.getWidth(), height = im.getHeight();

#else
    int width = 2048, height = 1489;
//    int width = 200, height = 200;
    printf("[no load]");
#endif
    printf("Loaded: %d x %d\n", width, height);

    //store image data in img_var(x, y, 0) and variance data in img_var(x, y, 1)
    Image<float> image(width, height);
    Image<float> variance(width, height);
    Image<uint16_t> mask(width, height);

#ifndef STANDALONE 
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
#endif
 
    int boundingBox = 2;
    Var x, y, c, i, yi, y0;


    Func polynomial ("polynomial");
    polynomial(x, y) = 0.1f + 0.001416015625f*x + 0.0f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.0f*x + 0.001416015625f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    Func kernel ("kernel"); 
    kernel(x, y, i) = exp(-(i*i)/(2*polynomial(x, y)*polynomial(x, y)));

    Func kernel1 ("kernel1"); 
    kernel1(x, y, i) = exp(-(i*i)/(2*polynomial1(x, y)*polynomial1(x, y)));


    Func in_bounded = BoundaryConditions::repeat_edge(image);

    //blur in the y direction
    Func blur_y ("blur_y");
    //does not compile with this order and vectorization
/*    Expr blur_image_helpY = 0.0f;
    Expr norm2 = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        blur_image_helpY += in_bounded(x, y + i) * kernel(x, y, i); 
        norm2 += kernel(x, y, i);
    }
    blur_image_helpY = blur_image_helpY/norm2;
    blur_y(x, y) = blur_image_helpY;
*/

    //does not compile with this order and vectorization
   blur_y(x, y) = (in_bounded(x, y + -2) * kernel1(x, y, -2)
                    + in_bounded(x, y + -1) * kernel1(x, y, -1)
                    + in_bounded(x, y + 0) * kernel1(x, y, 0)
                    + in_bounded(x, y + 1) * kernel1(x, y, 1)
                    + in_bounded(x, y + 2) * kernel1(x, y, 2))
                     /(kernel1(x, y, -1) + kernel1(x, y, 0) + kernel1(x, y, 1) 
                    + kernel1(x, y, 2) +  kernel1(x, y, -2));


    //compiles with this order
/*    blur_y(x, y) = (
                     in_bounded(x, y + -1) * kernel1(x, y, -1)
                    + in_bounded(x, y + -2) * kernel1(x, y, -2)
                    + in_bounded(x, y + 0) * kernel1(x, y, 0)
                    + in_bounded(x, y + 1) * kernel1(x, y, 1)
                    + in_bounded(x, y + 2) * kernel1(x, y, 2))
                     /(kernel1(x, y, -1) + kernel1(x, y, 0) + kernel1(x, y, 1) 
                    + kernel1(x, y, 2) +  kernel1(x, y, -2));
*/

    //blur in the x direction
    Func blur_x ("blur_x");
    //does not compile with this order
/*    Expr blur_image_helpX = 0.0f;
    Expr norm1 = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        blur_image_helpX += blur_y(x + i, y) * kernel(x, y, i); 
        norm1 += kernel(x, y, i);
    }
    blur_image_helpX = blur_image_helpX/norm1;
    blur_x(x, y) = blur_image_helpX;
*/
    blur_x(x, y) = (blur_y(x+ -1, y) * kernel(x, y, -1)
                    + blur_y(x + 0, y) * kernel(x, y, 0)
                    + blur_y(x + 1, y) * kernel(x, y, 1)
                    + blur_y(x + 2, y) * kernel(x, y, 2)
                    + blur_y(x + -2, y) * kernel(x, y, -2))
                     /(kernel(x, y, -1) + kernel(x, y, 0) + kernel(x, y, 1) 
                    + kernel(x, y, 2) +  kernel(x, y, -2));

    // Schedule
    //jit compilation hangs with vectorization:
    blur_x.compute_root().vectorize(x, 4).split(y, y0, yi, 16).parallel(y0);
    blur_y.compute_at(blur_x, y0).vectorize(x, 4);

    //jit compilation works without vectorization (with normal ordering of kernel terms in blur_y):
//    blur_x.compute_root().split(y, y0, yi, 16).parallel(y0);
//    blur_y.compute_at(blur_x, y0);

    polynomial.compute_root();



    blur_x.compile_to_lowered_stmt("blurFast.html", {image}, HTML);
    blur_x.print_loop_nest();
    // Benchmark the pipeline.
    Image<float> image_output(image.width(), image.height());
    //blur_x.realize(image_output);
    blur_x.compile_jit();
    cout << "jit compiled" << endl;

    Image<float> variance_output(variance.width(), variance.height());
//    blurVarianceY.realize(variance_output);

    Image<int32_t> mask_output(mask.width(), mask.height());
//    maskOutY.realize(mask_output);

    double average = 0;
    double min;
    double max;
    double imgTime;
    double varTime;
    double maskTime;
    int numberOfRuns = 5;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        blur_x.realize(image_output);
        double t2 = current_time();
//        blurVarianceY.realize(variance_output);
        double t3 = current_time();
//        maskOutY.realize(mask_output);
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
    "maskTime = " << maskTime << endl;

//    blur_mask.realize(mask_output);

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

    imOut.writeFits("./halideBlurFast.fits");
#endif

}


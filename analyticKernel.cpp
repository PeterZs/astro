// On os x:
//
//g++ -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -I ./include -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -L ./bin test.cpp -lafw -ldaf_base -ldaf_persistence -lboost_system -lHalide `libpng-config --cflags --ldflags` -o test
//
//DYLD_LIBRARY_PATH=DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/:./bin ./test
//
// On linux, you can compile and run it like so:
//g++ -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -I ./include -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -L ./bin test.cpp -lafw -ldaf_base -ldaf_persistence -lboost_system -lHalide `libpng-config --cflags --ldflags` -lpthread -ldl -o test
//
//LD_LIBRARY_PATH=DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/:./bin ./test
//
// On os x:
// g++ test.cpp -g -I ./include -L ./bin -lHalide `libpng-config --cflags --ldflags` -o test -std=c++11
// DYLD_LIBRARY_PATH=../bin ./test



// On os x:
// g++ analyticKernel.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o analyticKernel -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./analyticKernel



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

inline float cPolynomial(float x, float y){
    return (0.1f + 0.0f*x + 0.0019476158495634653f*y + 0.000001f*x*x + 0.000001f*x*y
    + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
    + 0.000000001f*y*y*y);
}

//check that this casts parameters from int to float and returns correct result
inline float cKernel(float x, float y, float i, float j, float maxX, float maxY){
    if(((x+i)>=maxX) || ((x+i)<0) || ((y+j)>=maxY) || ((y+j)<0))
        return 0;
    return ((exp(-((i*cos(cPolynomial(x, y)) +j*sin(cPolynomial(x, y)))
    *(i*cos(cPolynomial(x, y)) +j*sin(cPolynomial(x, y))))
    /(2*cPolynomial(x, y)*cPolynomial(x, y))) / (sqrtf(2*M_PI)*cPolynomial(x, y)))
    *(exp(-((j*cos(cPolynomial(x, y)) - i*sin(cPolynomial(x, y)))
    *(j*cos(cPolynomial(x, y)) - i*sin(cPolynomial(x, y))))
    /(2*cPolynomial(x, y)*cPolynomial(x, y))) / (sqrtf(2*M_PI)*cPolynomial(x, y))));
}

//return 16 bits all set to 1 when
//the kernel is not 0 and 16 bits all set to 0 when the kernel is 0
/*inline uint16_t maskHelper(int x, int y, int i, int j){
    return (~(((uint16_t)(cKernel(x, y, i, j) != 0.0f)) - 1));
}
*/

int main(int argc, char *argv[]) {
	auto im = afwImage::MaskedImage<float>("../calexp-004207-g3-0123.fits");
    printf("Loaded: %d x %d\n", im.getWidth(), im.getHeight());

    //store image data in img_var(x, y, 0) and variance data in img_var(x, y, 1)
    Image<float> img_var(im.getWidth(), im.getHeight(), 2);
    Image<uint16_t> mask(im.getWidth(), im.getHeight());

    //Read image in
    for (int y = 0; y < im.getHeight(); y++) {
    	afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = im.x_at(0, y);
	   	for (int x = 0; x < im.getWidth(); x++){
       		img_var(x, y, 0) = (*inPtr).image();
      		img_var(x, y, 1) = (*inPtr).variance();
     		mask(x, y) = (*inPtr).mask();
     		inPtr++;
        }
    }


    int boundingBox = 2; 
    Var x, y, i, j, i_v, y0, yi;

    //compute output image and variance
    Func polynomial;
    polynomial(x, y) = 0.1f + 0.0f*x + 0.0019476158495634653f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //The kernel is spatially variant for different positions (x, y) in the image
    //kernel(x, y, i, j) is the value of the kernel at position (x, y) in the image
    //and position (i, j) inside the kernel
    Func kernel;
    kernel(x, y, i, j) = (exp(-((i*cos(polynomial(x, y)) +j*sin(polynomial(x, y)))
                    *(i*cos(polynomial(x, y)) +j*sin(polynomial(x, y))))
                    /(2*polynomial(x, y)*polynomial(x, y))) / (sqrtf(2*M_PI)*polynomial(x, y)))
                    *(exp(-((j*cos(polynomial(x, y)) - i*sin(polynomial(x, y)))
                    *(j*cos(polynomial(x, y)) - i*sin(polynomial(x, y))))
                    /(2*polynomial(x, y)*polynomial(x, y))) / (sqrtf(2*M_PI)*polynomial(x, y)));

    //calculate the spatially dependent normalization
    Func norm;
    Expr norm_help = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            norm_help += (exp(-((i*cos(polynomial(x, y)) +j*sin(polynomial(x, y)))
                    *(i*cos(polynomial(x, y)) +j*sin(polynomial(x, y))))
                    /(2*polynomial(x, y)*polynomial(x, y))) / (sqrtf(2*M_PI)*polynomial(x, y)))
                    *(exp(-((j*cos(polynomial(x, y)) - i*sin(polynomial(x, y)))
                    *(j*cos(polynomial(x, y)) - i*sin(polynomial(x, y))))
                    /(2*polynomial(x, y)*polynomial(x, y))) / (sqrtf(2*M_PI)*polynomial(x, y)));
        }
    }
    norm(x, y) = norm_help;

    Func img_var_bounded = BoundaryConditions::repeat_edge(img_var);

    Func blur;
    blur(x, y, i_v) = 0.0f;
    //compute Image output
    Expr blur_image_help = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            blur_image_help += img_var_bounded(x + i, y + j, 0) * kernel(x, y, i, j); 
        }
    }
    blur_image_help = blur_image_help/norm(x, y);
    blur(x, y, 0) = blur_image_help;


    //compute Variance output
    Expr blur_variance_help = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            blur_variance_help += img_var_bounded(x + i, y + j, 1) * kernel(x, y, i, j) * kernel(x, y, i, j); 
        }
    }
    blur_variance_help = blur_variance_help/(norm(x,y)*norm(x,y));
    blur(x, y, 1) = blur_variance_help;

    //Schedule
  //  blur.reorder(i_v, x, y);
    // Split the y coordinate of the consumer into strips of 16 scanlines:
    blur.split(y, y0, yi, 40);
    // Compute the strips using a thread pool and a task queue.
    blur.parallel(y0);
    // Vectorize across x by a factor of four.
    blur.vectorize(x, 8);

//    polynomial.compute_at(blur, x).vectorize(x, 8);


    // Print out pseudocode for the pipeline.
    blur.compile_to_lowered_stmt("blur.html", {img_var}, HTML);


    //testing
/*    bool t1 = false;
    uint16_t test1 = (uint16_t) t1;
    test1 = ~(test1-1);
    bitset<16> test1b(test1);

    bool t2 = true;
    uint16_t test2 = (uint16_t) t2;
    test2 = ~(test2-1);
    bitset<16> test2b(test2);

    cout << "test1 = " << test1 << ", test2 = " << test2 << endl;
*/
    //done testing




    //compute output mask
    //~(((uint16_t)(kernel(x, y, i, j) != 0.0f)) - 1) should produce 16 bits all set to 1 when
    //the kernel is not 0 and 16 bits all set to 0 when the kernel is 0
    //the new mask is the bitwise OR of each mask whose kernel value is nonzero
   //Current Implementation
/*    Expr blur_mask_help = 0;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
                blur_mask_help = ((~(((uint16_t)(kernel(x, y, i, j) != 0.0f)) - 1) & mask_bounded(x + i, y + j)) | blur_mask_help); 
        }
    }
    Func blur_mask;
    blur_mask(x, y) = blur_mask_help;
*/
    //temp SUPER slow, don't use update, try image instead
/*    Func blur_mask;
    blur_mask(x, y) = 0;
    for(int y1 = 0; y1 < mask.height(); y1++){
        for(int x1 = 0; x1 < mask.width(); x1++){
            for(int i = -boundingBox; i <= boundingBox; i++){
                for(int j = -boundingBox; j <= boundingBox; j++){
                    blur_mask(x1, y1) = blur_mask(x1, y1) | (maskHelper(x1, y1, i, j) & mask_bounded(x1 + i, y1 + j));
                }
            }
        }
        cout << "finished row:" << y1 << endl;
    }
*/    
    //done temp
    //Not halide
 /*   for(int x1 = 0; x1 < mask.width(); x1++){
        for(int y1 = 0; y1 < mask.height(); y1++){
            for(int i = -boundingBox; i <= boundingBox; i++){
                for(int j = -boundingBox; j <= boundingBox; j++){

                }
            }
        }
    }
*/

    //Compute mask
    Func mask_bounded = BoundaryConditions::repeat_edge(mask);

    Func maskOut;

    Expr maskOutHelp = 0;

    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
//            maskOutHelp = select(kernel(x, y, i, j) == 0.0f, maskOutHelp, maskOutHelp | mask_bounded(x + i, y + j));
            maskOutHelp = maskOutHelp | mask_bounded(x + i, y + j);

/*                if((x1 < 10) && (y1 < 10)){
                    cout << "(x, y) = (" << x1 << ", " << y1 << ")";
                    cout << "(i, j) = (" << i << ", " << j << ")";
                    cout << " kernel = " << cKernel(x1, y1, i, j, mask.width(), mask.height()) << ", ";
                    cout << " mask_output = " << mask_output(x1, y1) << ", ";
                    cout << " T/F?: " << (cKernel(x1, y1, i, j, mask.width(), mask.height()) != 0.0f) << endl;
                }
*/      
        }
    }

    

    maskOut(x, y) = maskOutHelp;

    // Split the y coordinate of the consumer into strips of 16 scanlines:
    maskOut.split(y, y0, yi, 30);
    // Compute the strips using a thread pool and a task queue.
    maskOut.parallel(y0);
    // Vectorize across x by a factor of four.
    maskOut.vectorize(x, 8);

    
    // Benchmark the pipeline.
    Image<float> img_var_output(img_var.width(), img_var.height(), 2);
    blur.realize(img_var_output);

    Image<int32_t> mask_output(mask.width(), mask.height());
    maskOut.realize(mask_output);

    double average = 0;
    double min;
    double max;
    int numberOfRuns = 5;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        blur.realize(img_var_output);
        maskOut.realize(mask_output);
        double t2 = current_time();
        double curTime = (t2-t1);
        average += curTime;
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
    average = average/numberOfRuns;
    std::cout << "Average Time: " << average << ", Min = " <<
    min << ", Max = " << max << ", with " << numberOfRuns <<
    " runs" << '\n';

//    blur_mask.realize(mask_output);

    //write image out
    auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
    for (int y = 0; y < imOut.getHeight(); y++) {
    	afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);

        for (int x = 0; x < imOut.getWidth(); x++){
        	afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
            curPixel(img_var_output(x, y, 0), mask_output(x, y), img_var_output(x, y, 1));
        	(*inPtr) = curPixel;
        	inPtr++;

        }
    }

	imOut.writeFits("./halideAnalyticKernel.fits");
}


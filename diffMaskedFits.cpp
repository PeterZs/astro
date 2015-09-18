// On os x:
// g++ diffMaskedFits.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o diffMaskedFits -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./diffMaskedFits
// On linux:
//g++ diffMaskedFits.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o diffMaskedFits -std=c++11
//
//LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./diffMaskedFits

//for the purpose of creating output image of (image diff)/trueVariance, enter reference image as first argument (its variance will be used)

#include "lsst/afw/image.h"
#include <stdio.h>

#include "Halide.h"
using namespace std;
using namespace Halide;

using Halide::Image;

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;


int main(int argc, char *argv[]) {
    if(argc != 3){
        cout << "error, must supply two .fits files as arguments" << endl;
        return 0;
    }
	auto im1 = afwImage::MaskedImage<float>(argv[1]);
    printf("Loaded: %d x %d\n", im1.getWidth(), im1.getHeight());
    Image<float> image1(im1.getWidth(), im1.getHeight());
    Image<float> variance1(im1.getWidth(), im1.getHeight());
    Image<uint16_t> mask1(im1.getWidth(), im1.getHeight());

    auto im2 = afwImage::MaskedImage<float>(argv[2]);
    printf("Loaded: %d x %d\n", im2.getWidth(), im2.getHeight());
    Image<float> image2(im2.getWidth(), im2.getHeight());
    Image<float> variance2(im2.getWidth(), im2.getHeight());
    Image<uint16_t> mask2(im2.getWidth(), im2.getHeight());


    //Read image1 in
    for (int y = 0; y < im1.getHeight(); y++) {
    	afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>::x_iterator inPtr = im1.x_at(0, y);
	   	for (int x = 0; x < im1.getWidth(); x++){
       		image1(x, y) = (*inPtr).image();
      		variance1(x, y) = (*inPtr).variance();
     		mask1(x, y) = (*inPtr).mask();
     		inPtr++;
        }
    }

    //Read image2 in
    for (int y = 0; y < im2.getHeight(); y++) {
        afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>::x_iterator inPtr = im2.x_at(0, y);
        for (int x = 0; x < im2.getWidth(); x++){
            image2(x, y) = (*inPtr).image();
            variance2(x, y) = (*inPtr).variance();
            mask2(x, y) = (*inPtr).mask();
            inPtr++;
        }
    }

    if((im1.getWidth() != im2.getWidth()) || (im2.getHeight() != im2.getHeight())){
        cout << "image sizes do not match" << endl;
        return 0;
    }

    //compute the maximum differences of the input images' image, mask, and variance values
    //write an output image whose values are image1-image2 
    double maxImageDiff = 0;
    double maxImageDiffPercent = 0;
    double valueAtMaxImageDiff = 0;
    double varianceAtMaxImgDiff = 0;
    int imageX = -1;
    int imageY = -1;

    double maxVarianceDiff = 0;
    double maxVarianceDiffPercent = 0;
    double valueAtMaxVarianceDiff1 = 0;
    double valueAtMaxVarianceDiff2 = 0;
    int varianceX = -1;
    int varianceY = -1;

    int maxMaskDiff = 0;
    int maskX = -1;
    int maskY = -1;
    int mask1Val = -1;
    int mask2Val = -1;

    float maxImage1 = 0;
    float maxVariance1 = 0;
    int maxMask1 = 0;

    float maxImage2 = 0;
    float maxVariance2 = 0;
    int maxMask2 = 0;


    int kernelSize = 10;
    auto imDiff = afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>(im1.getWidth(), im1.getHeight());
    for (int y = kernelSize; y < imDiff.getHeight()-kernelSize; y++) {
    	afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>::x_iterator inPtr = imDiff.x_at(0, y);

        for (int x = kernelSize; x < imDiff.getWidth()-kernelSize; x++){
//            afwImage::pixel::SinglePixel<float, afwImage::MaskPixel, afwImage::VariancePixel> curPixel(image1(x, y) - image2(x, y), variance1(x, y) - variance2(x, y), mask1(x, y) - mask2(x, y));
            afwImage::pixel::SinglePixel<float, afwImage::MaskPixel, afwImage::VariancePixel> curPixel((image1(x, y) - image2(x, y))/variance1(x, y), variance1(x, y) - variance2(x, y), mask1(x, y) - mask2(x, y));
        	(*inPtr) = curPixel;
        	inPtr++;
//            if(abs(image1(x, y) - image2(x, y))/min(abs(image1(x, y)), abs(image2(x, y))) > maxImageDiffPercent){
            if(abs(image1(x, y) - image2(x, y)) > maxImageDiff){
                maxImageDiffPercent = abs(image1(x, y) - image2(x, y))/min(abs(image1(x, y)), abs(image2(x, y)));
                maxImageDiff = abs(image1(x, y) - image2(x, y));
                valueAtMaxImageDiff = min(abs(image1(x, y)), abs(image2(x, y)));
                varianceAtMaxImgDiff = min(abs(variance1(x, y)), abs(variance2(x, y)));
                imageX = x;
                imageY = y;
            }
            if(abs(variance1(x, y) - variance2(x, y))/min(abs(variance1(x, y)), abs(variance2(x, y))) > maxVarianceDiff){
                maxVarianceDiffPercent = abs(variance1(x, y) - variance2(x, y))/min(abs(variance1(x, y)), abs(variance2(x, y)));
                maxVarianceDiff = abs(variance1(x, y) - variance2(x, y));
                valueAtMaxVarianceDiff1 = variance1(x, y);
                valueAtMaxVarianceDiff2 = variance2(x, y);
                varianceX = x;
                varianceY = y;
            }
            if(abs(mask1(x, y) - mask2(x, y)) > maxMaskDiff){
                maxMaskDiff = abs(mask1(x, y) - mask2(x, y));
                maskX = x;
                maskY = y;
                mask1Val = mask1(x, y);
                mask2Val = mask2(x, y);
            }

            if(abs(image1(x, y)) > maxImage1)
                maxImage1 = abs(image1(x, y));
            if(abs(variance1(x, y)) > maxVariance1)
                maxVariance1 = abs(variance1(x, y));
            if(mask1(x, y) > maxMask1)
                maxMask1 = mask1(x, y);

            if(abs(image2(x, y)) > maxImage2)
                maxImage2 = abs(image2(x, y));
            if(abs(variance2(x, y)) > maxVariance2)
                maxVariance2 = abs(variance2(x, y));
            if(mask2(x, y) > maxMask2)
                maxMask2 = mask2(x, y);
        }
    }

    cout << "Max (image difference)/(min img value) = " << maxImageDiffPercent << ",  Max image difference = " << maxImageDiff
    << ", value at max image difference = " << valueAtMaxImageDiff << " at position: (" << imageX << ", " << imageY
    << "), " << "variance at max difference = " << varianceAtMaxImgDiff << endl;

    cout << "Max (variance difference)/(min var value) = " << maxVarianceDiffPercent << ",  Max Variance difference = " << maxVarianceDiff
    << ", variance1 at max Variance difference = " << valueAtMaxVarianceDiff1 << ", variance2 at max Variance difference = "
    << valueAtMaxVarianceDiff2 <<" at position: (" << varianceX << ", ";
    cout << varianceY << ")" << endl;

    cout << "Max mask difference = " << maxMaskDiff << " at position: (" << maskX << ", " << maskY << ")" << 
    ", img1 mask = " << mask1Val << ", img2 mask = " << mask2Val << endl;

    cout << "Max image1 = " << maxImage1 << ", max variance1 = " << maxVariance1 << 
    ", max mask1 = " << maxMask1 << endl;

    cout << "Max image2 = " << maxImage2 << ", max variance2 = " << maxVariance2 << 
    ", max mask2 = " << maxMask2 << endl;

//	imDiff.writeFits("./imageDifferenceDivVariance.fits");
}


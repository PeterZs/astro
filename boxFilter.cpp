// Download a Halide distribution from halide-lang.org and untar it in
// the current directory. Then you should be able to compile this
// file with:
//
//
//Using FITS images on os x:
// g++ Gaussian2D.cpp -g -I ../include -L ../bin -lHalide -lcfitsio `libpng-config --cflags --ldflags` -o Gaussian2D -std=c++11
// DYLD_LIBRARY_PATH=../bin ./Gaussian2D
//
//Using FITS images on linux:
// g++ Gaussian2D.cpp -g -I ../include -L ../bin -lHalide -lpthread -ldl -lcfitsio `libpng-config --cflags --ldflags` -o Gaussian2D -std=c++11
// LD_LIBRARY_PATH=../bin ./Gaussian2D
//
//
// You'll also need a multi-megapixel png image to run this on. Name
// it input.png and put it in this directory.

// Include the Halide language


// On os x:
// g++ boxFilter.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o boxFilter -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./boxFilter

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

int main(int argc, char **argv) {
	auto im = afwImage::MaskedImage<float>("../calexp-004207-g3-0123.fits");
    printf("Loaded: %d x %d\n", im.getWidth(), im.getHeight());

    //store image data in img_var(x, y)
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
	
	int boundingBox = 2;
    
    Var x, y, c, y0, yi, i, j;

	//Kernel #1
    Func totalKernel;
    float sigmaX1 = 200.0f;
    float sigmaY1 = 200.0f;
	float theta1 = 0.0f; //rotation of sigmaX axis
    totalKernel(i, j) = (exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
					/(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
					*(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
					/(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1));




	float norm = 0.0f;
	float curKernelPos;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			curKernelPos = 
					((exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
					/(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
					*(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
					/(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1)))
					;
			norm += curKernelPos;

			//debugging
			cout << curKernelPos << "\t";
		}
		//debugging
		cout << endl;
	}
	cout << "Norm = " << norm << endl;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			curKernelPos = 
					((exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
					/(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
					*(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
					/(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1)))
					;
			//debugging
			cout << curKernelPos/norm << "\t";
		}
		//debugging
		cout << endl;
	}

    Func image_bounded = BoundaryConditions::repeat_edge(image);
    Func variance_bounded = BoundaryConditions::repeat_edge(variance);

    //compute new Image values
	Func blur_image;
	Expr blur_image_help = 0.0f;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			blur_image_help += image_bounded(x + i, y + j) * totalKernel(i, j); 
		}
	}
	blur_image_help = blur_image_help/norm;

	blur_image(x, y) = blur_image_help;

	//compute new Variance values 

	Func blur_variance;
	Expr blur_variance_help = 0.0f;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			blur_variance_help += variance_bounded(x + i, y + j)*totalKernel(i, j)/norm * totalKernel(i, j)/norm; 
//			blur_variance_help += variance_bounded(x + i, y + j)/(25*25); 
//			blur_variance_help += variance_bounded(x + i, y + j) * totalKernel(i, j)/norm; 


		}
	}

	blur_variance(x, y) = blur_variance_help;

//	Func blur_variance;
//	Expr blur_variance_help = 0.0f;
//	for(int i = -boundingBox; i <= boundingBox; i++){
//		for(int j = -boundingBox; j <= boundingBox; j++){
//			blur_variance_help += totalKernel(i, j)/norm * 
//			(variance_bounded(x + i, y + j) + image_bounded(x + i, y + j)*image_bounded(x + i, y + j)); 
//		}
//	}
//	blur_variance_help = blur_variance_help - blur_image_help*blur_image_help;
//
//	blur_variance(x, y) = blur_variance_help;
//*****************************************************************************
//	Func blur_variance;
//	Expr curVarianceMean = 0.0f;
//	for(int i = -boundingBox; i <= boundingBox; i++){
//		for(int j = -boundingBox; j <= boundingBox; j++){
//			curVarianceMean += totalKernel(i, j)/norm * variance_bounded(x + i, y + j); 
//		}
//	}
//	Expr varianceOfVariance = 0.0f;
//	for(int i = -boundingBox; i <= boundingBox; i++){
//		for(int j = -boundingBox; j <= boundingBox; j++){
//			varianceOfVariance += totalKernel(i, j)/norm * (variance_bounded(x + i, y + j) - curVarianceMean)
//			* (variance_bounded(x + i, y + j) - curVarianceMean); 
//		}
//	}
//	blur_variance(x, y) = varianceOfVariance;

	//end variance

    // Split the y coordinate of the consumer into strips of 16 scanlines:
    blur_image.split(y, y0, yi, 30);
    blur_variance.split(y, y0, yi, 30);

    // Compute the strips using a thread pool and a task queue.
    blur_image.parallel(y0);
    blur_variance.parallel(y0);

    // Vectorize across x by a factor of four.
    blur_image.vectorize(x, 8);
    blur_variance.vectorize(x, 8);

//	for(int i = 0; i < boundingBox*boundingBox; i++)
//		blur.update(i).vectorize(x, 8);


    // Print out pseudocode for the pipeline.
//




//    blur.compile_to_lowered_stmt("blur_image.html", {image}, HTML);
    
    // Benchmark the pipeline.
    Image<float> image_output(image.width(), image.height());
    Image<float> variance_output(image.width(), image.height());

    blur_image.realize(image_output);
    blur_variance.realize(variance_output);



    //Compute mask
    Func mask_bounded = BoundaryConditions::repeat_edge(mask);

    Func maskOut;

    Expr maskOutHelp = 0;
    for(int i = -boundingBox; i <= boundingBox; i++)
    	for(int j = -boundingBox; j <= boundingBox; j++){
    		//bitwise OR in following line seems to cast result to int32, so mask_output must be typed int32, is this a problem?
    		maskOutHelp = select(totalKernel(i, j) == 0.0f, maskOutHelp, maskOutHelp | mask_bounded(x + i, y + j));


    		//maskOutHelp = select(totalKernel(i, j) == 0.0f, maskOutHelp, maskOutHelp + mask_bounded(x + i, y + j));

    	}

    maskOut(x, y) = maskOutHelp;

    // Split the y coordinate of the consumer into strips of 16 scanlines:
    maskOut.split(y, y0, yi, 30);
    // Compute the strips using a thread pool and a task queue.
    maskOut.parallel(y0);
    // Vectorize across x by a factor of four.
    maskOut.vectorize(x, 8);

    //
    Image<int32_t> mask_output(mask.width(), mask.height());

    maskOut.realize(mask_output);


    double t1 = current_time();
    maskOut.realize(mask_output);
    double t2 = current_time();
    cout << "mask computation took " << t2-t1 << endl;

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

	imOut.writeFits("./halideBox.fits");

/*	
	double average = 0;
	double min;
	double max;
	int numberOfRuns = 20;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        blur.realize(output);
        double t2 = current_time();
//        std::cout << "Time: " << (t2 - t1) << '\n';
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
   */


  
    return 0;
}





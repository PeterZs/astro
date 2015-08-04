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
// g++ spatiallyInvariantGaussianKernel.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o spatiallyInvariantGaussianKernel -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./spatiallyInvariantGaussianKernel

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
    
    Var x, y, c, y0, yi, i, j;

	//Kernel #1
    Func kernel1;
    float sigmaX1 = 2.0f;
    float sigmaY1 = 2.0f;
	float theta1 = 0.0f; //rotation of sigmaX axis
    kernel1(i, j) = (exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
					/(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
					*(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
					/(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1));



	//Kernel #2
    Func kernel2;
    float sigmaX2 = 0.5f;
    float sigmaY2 = 4.0f;
	float theta2 = 0.0f; //rotation of sigmaX axis
    kernel2(i, j) = (exp(-((i*cos(theta2) +j*sin(theta2))*(i*cos(theta2) +j*sin(theta2)))
					/(2*sigmaX2*sigmaX2)) / (sqrtf(2*M_PI)*sigmaX2))
					*(exp(-((j*cos(theta2) - i*sin(theta2))*(j*cos(theta2) - i*sin(theta2)))
					/(2*sigmaY2*sigmaY2)) / (sqrtf(2*M_PI)*sigmaY2));

	//Kernel #3
    Func kernel3;
    float sigmaX3 = 0.5f;
    float sigmaY3 = 4.0f;
	float theta3 = 3.14159f/4; //rotation of sigmaX axis
    kernel3(i, j) = (exp(-((i*cos(theta3) +j*sin(theta3))*(i*cos(theta3) +j*sin(theta3)))
					/(2*sigmaX3*sigmaX3)) / (sqrtf(2*M_PI)*sigmaX3))
					*(exp(-((j*cos(theta3) - i*sin(theta3))*(j*cos(theta3) - i*sin(theta3)))
					/(2*sigmaY3*sigmaY3)) / (sqrtf(2*M_PI)*sigmaY3));
	//Kernel #4
    Func kernel4;
    float sigmaX4 = 0.5f;
    float sigmaY4 = 4.0f;
	float theta4 = 3.14159f/2; //rotation of sigmaX axis
    kernel4(i, j) = (exp(-((i*cos(theta4) +j*sin(theta4))*(i*cos(theta4) +j*sin(theta4)))
					/(2*sigmaX4*sigmaX4)) / (sqrtf(2*M_PI)*sigmaX4))
					*(exp(-((j*cos(theta4) - i*sin(theta4))*(j*cos(theta4) - i*sin(theta4)))
					/(2*sigmaY4*sigmaY4)) / (sqrtf(2*M_PI)*sigmaY4));


	//Kernel #5
    Func kernel5;
    float sigmaX5 = 4.0f;
    float sigmaY5 = 4.0f;
	float theta5 = 0.0; //rotation of sigmaX axis
    kernel5(i, j) = (exp(-((i*cos(theta5) +j*sin(theta5))*(i*cos(theta5) +j*sin(theta5)))
					/(2*sigmaX5*sigmaX5)) / (sqrtf(2*M_PI)*sigmaX5))
					*(exp(-((j*cos(theta5) - i*sin(theta5))*(j*cos(theta5) - i*sin(theta5)))
					/(2*sigmaY5*sigmaY5)) / (sqrtf(2*M_PI)*sigmaY5));

	Func totalKernel;
	totalKernel(i, j) = kernel1(i, j) + kernel2(i, j) + kernel3(i, j) + kernel4(i, j) + kernel5(i, j);

	float norm = 0.0f;
	float curKernelPos;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			curKernelPos = 
					((exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
					/(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
					*(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
					/(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1))) +

					((exp(-((i*cos(theta2) +j*sin(theta2))*(i*cos(theta2) +j*sin(theta2)))
					/(2*sigmaX2*sigmaX2)) / (sqrtf(2*M_PI)*sigmaX2))
					*(exp(-((j*cos(theta2) - i*sin(theta2))*(j*cos(theta2) - i*sin(theta2)))
					/(2*sigmaY2*sigmaY2)) / (sqrtf(2*M_PI)*sigmaY2))) + 

					((exp(-((i*cos(theta3) +j*sin(theta3))*(i*cos(theta3) +j*sin(theta3)))
					/(2*sigmaX3*sigmaX3)) / (sqrtf(2*M_PI)*sigmaX3))
					*(exp(-((j*cos(theta3) - i*sin(theta3))*(j*cos(theta3) - i*sin(theta3)))
					/(2*sigmaY3*sigmaY3)) / (sqrtf(2*M_PI)*sigmaY3))) + 

					((exp(-((i*cos(theta4) +j*sin(theta4))*(i*cos(theta4) +j*sin(theta4)))
					/(2*sigmaX4*sigmaX4)) / (sqrtf(2*M_PI)*sigmaX4))
					*(exp(-((j*cos(theta4) - i*sin(theta4))*(j*cos(theta4) - i*sin(theta4)))
					/(2*sigmaY4*sigmaY4)) / (sqrtf(2*M_PI)*sigmaY4))) +

					((exp(-((i*cos(theta5) +j*sin(theta5))*(i*cos(theta5) +j*sin(theta5)))
					/(2*sigmaX5*sigmaX5)) / (sqrtf(2*M_PI)*sigmaX5))
					*(exp(-((j*cos(theta5) - i*sin(theta5))*(j*cos(theta5) - i*sin(theta5)))
					/(2*sigmaY5*sigmaY5)) / (sqrtf(2*M_PI)*sigmaY5)))
					;
			norm += curKernelPos;

			//debugging
			cout << curKernelPos << "\t";
		}
		//debugging
		cout << endl;
	}
	cout << "Norm = " << norm << endl;

    Func in_bounded = BoundaryConditions::repeat_edge(img_var);

	Func blur;
	blur(x, y, c) = 0.0f;
    //compute Image output
	Expr blur_image_help = 0.0f;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			blur_image_help += in_bounded(x + i, y + j, 0) * totalKernel(i, j); 
		}
	}
	blur_image_help = blur_image_help/norm;
	blur(x, y, 0) = blur_image_help;


	//compute Variance output
	Expr blur_variance_help = 0.0f;
	for(int i = -boundingBox; i <= boundingBox; i++){
		for(int j = -boundingBox; j <= boundingBox; j++){
			blur_variance_help += in_bounded(x + i, y + j, 1) * totalKernel(i, j) * totalKernel(i, j); 
		}
	}
	blur_variance_help = blur_variance_help/(norm*norm);
	blur(x, y, 1) = blur_variance_help;

	blur.reorder(c, x, y);
    // Split the y coordinate of the consumer into strips of 16 scanlines:
    blur.split(y, y0, yi, 30);
    // Compute the strips using a thread pool and a task queue.
    blur.parallel(y0);
    // Vectorize across x by a factor of four.
    blur.vectorize(x, 8);

//	for(int i = 0; i < boundingBox*boundingBox; i++)
//		blur.update(i).vectorize(x, 8);


    // Print out pseudocode for the pipeline.
    blur.compile_to_lowered_stmt("blur.html", {img_var}, HTML);


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

	imOut.writeFits("./halideSpatiallyInvariantGaussian.fits");

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





#ifndef __LINEAR_COMBINATION_KERNEL_1__
#define __LINEAR_COMBINATION_KERNEL_1__


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


#define BOUNDING_BOX 2 //the kernel has dimensions (BOUNDING_BOX*2 + 1) x (BOUNDING_BOX*2 + 1)
#define NUMBER_OF_RUNS 5 //number of runs when performance testing

// from lesson_12_using_the_gpu.cpp
// We're going to want to schedule a pipeline in several ways, so we
// define the pipeline in a class so that we can recreate it several
// times with different schedules.

// Define some Vars to use.
Var x, y, i, j, i_v, y_0, yi;

//CURRENTLY ONLY IMAGE PLANE IS IMPLEMENTED
//this program uses tuples
//each basis kernel is convolved with each input plane (using tuples) and then
//the 5 output planes are combined using the weights of the spatially varying
//polynomials
class convolveKernelsSeparatelyThenCombinePipeline {
public:
    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, kernel1, kernel2,
    kernel3, kernel4, kernel5, image_bounded, blurImage1, blurImage2, blurImage3,
    blurImage4, blurImage5, variance_bounded, mask_bounded, combined_output, imageOut,
    varianceOut, maskOut;

    Image<float> image;
    Image<float> variance;
    Image<uint16_t> mask;

    Buffer image_gpu_output;
    Buffer variance_gpu_output;
    Buffer mask_gpu_output;

    convolveKernelsSeparatelyThenCombinePipeline(Image<float> image_, Image<float> variance_,
        Image<uint16_t> mask_, bool useTuples_);

    void debug();

    void schedule_for_cpu();

    void schedule_for_gpu();

    void test_performance_cpu();

    void test_performance_gpu();

private:
    bool useTuples;

};

//this program uses tuples
//the total kernel is computed as a spatially varying linear combination of each basis
//kernel and convolved once with each input plane (using a tuple)
class convolveOneSpatiallyVaryingKernelPipeline {
public:
    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, kernel1, kernel2,
    kernel3, kernel4, kernel5, blurImage1, blurImage2, blurImage3, blurImage4, blurImage5,
    image_bounded, variance_bounded, mask_bounded, combined_output, imageOut, varianceOut,
    maskOut;

    Image<float> image;
    Image<float> variance;
    Image<uint16_t> mask;

    Buffer image_gpu_output;
    Buffer variance_gpu_output;
    Buffer mask_gpu_output;

    //will use tuples if useTuples==true, will not use tuples if useTuples==false
    convolveOneSpatiallyVaryingKernelPipeline(Image<float> image_, Image<float> variance_,
        Image<uint16_t> mask_, bool useTuples_);

    void schedule_for_cpu();

	void schedule_for_gpu();

    void test_performance_cpu();

    void test_performance_gpu();

    //check whether the image planes match
    //implement more if desired later
    void test_correctness(Image<float> reference_output);

    void debug();

private:
    //will use tuples if useTuples==true, will not use tuples if useTuples==false
    bool useTuples;
};


#endif //__LINEAR_COMBINATION_KERNEL_H__
#ifndef __LINEAR_COMBINATION_KERNEL_1__
#define __LINEAR_COMBINATION_KERNEL_1__


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
// Define some Vars to use.
Var x, y, i, j, i_v, y_0, yi;


// from lesson_12_using_the_gpu.cpp
// We're going to want to schedule a pipeline in several ways, so we
// define the pipeline in a class so that we can recreate it several
// times with different schedules.
class convolveKernelsSeparatelyThenCombinePipeline {
public:
    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, kernel1, kernel2,
    kernel3, kernel4, kernel5, image_bounded, blurImage1, blurImage2, blurImage3,
    blurImage4, blurImage5, variance_bounded, mask_bounded, combined_output;

    Image<float> image;
    Image<float> variance;
    Image<uint16_t> mask;

    Image<float> image_output;
    Image<float> variance_output;
    Image<uint16_t> mask_output;

    convolveKernelsSeparatelyThenCombinePipeline(Image<float> image_, Image<float> variance_,
        Image<uint16_t> mask_);

    void debug();

    void schedule_for_cpu();

    void schedule_for_gpu();

    void test_performance_cpu();

    void test_performance_gpu();

};


class convolveOneSpatiallyVaryingKernelPipeline {
public:
    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, kernel1, kernel2,
    kernel3, kernel4, kernel5, image_bounded, blurImage1, blurImage2, blurImage3,
    blurImage4, blurImage5, variance_bounded, mask_bounded, combined_output;

    Image<float> image;
    Image<float> variance;
    Image<uint16_t> mask;

    Image<float> image_output;
    Image<float> variance_output;
    Image<uint16_t> mask_output;

    convolveOneSpatiallyVaryingKernelPipeline(Image<float> image_, Image<float> variance_,
        Image<uint16_t> mask_);

    void schedule_for_cpu();

	void schedule_for_gpu();

    void test_performance_cpu();

    void test_performance_gpu();


    void debug();
};


#endif //__META_PATH_OUTLIER_H__
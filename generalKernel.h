#ifndef __GENERAL_KERNEL__
#define __GENERAL_KERNEL__


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
#include <math.h>       
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

//the total kernel is computed as a (possibly spatially varying) linear combination of each basis
//kernel (if there are multiple basis kernels) and convolved once with each input plane.
//Using a tuple to combine the computation of the three input planes is optional.  Tuples are
//generally faster on the CPU but slower on the GPU.
class convolveOncePipeline {
public:
    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, kernel1, kernel2,
    kernel3, kernel4, kernel5, blurImage1, blurImage2, blurImage3, blurImage4, blurImage5,
    image_bounded, variance_bounded, mask_bounded, combined_output, imageOut, varianceOut,
    maskOut;

    Image<float> image;
    Image<float> variance;
    Image<uint16_t> mask;

    Image<float> image_output;
    Image<float> variance_output;
    Image<uint16_t> mask_output;

    //will use tuples if useTuples==true, will not use tuples if useTuples==false
    convolveOncePipeline(Image<float> image_, Image<float> variance_,
        Image<uint16_t> mask_, bool useTuples_);

    virtual void schedule_for_cpu();

	virtual void schedule_for_gpu();

    virtual void test_performance_cpu();

    virtual void test_performance_gpu();

    //check whether the image planes match
    //implement more if desired later
    void test_correctness(Image<float> reference_output);

    void debug();

private:
    //will use tuples if useTuples==true, will not use tuples if useTuples==false
    bool useTuples;
};


//represents polynomial:
//polynomial(x, y) = a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y 
//                    uvv*x*y*y + vvv*y*y*y
class polynomial{
public:
    polynomial(float a_, float u_, float v_, float uu_, float uv_, float vv_, float uuu_,
        float uuv_, float uvv_, float vvv_): a(a_),  u(u_),  v(v_),  uu(uu_),  
        uv(uv_), vv(vv_),  uuu(uuu_), uuv(uuv_),  uvv(uvv_),  vvv(vvv_) {}

    Halide::Expr operator()(Halide::Var x, Halide::Var y){
        return (a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y +
                    uvv*x*y*y + vvv*y*y*y);
    }


private:
    float a, u, v, uu, uv, vv, uuu, uuv, uvv, vvv;

};


//Represents a 2 dimensional gaussian with standard deviations sigmaI and sigmaJ in its
//i and j dimensions respectively.  The guassian's i, j coordinate system is rotated by theta 
//radians with respect to the image's x, y coordinate system.
//The guassian is not normalized.
class gaussian2D{
public:
    gaussian2D(float sigmaI_, float sigmaJ_, float theta_)
        : sigmaI(sigmaI_), sigmaJ(sigmaJ_), theta(theta_){}

    gaussian2D(){} //empty constructor for use by child

    Halide::Expr operator()(Halide::Var x, Halide::Var y){
        return (exp(-((x*cos(theta) + y*sin(theta)) * (x*cos(theta) + y*sin(theta)))
                    /(2*sigmaI*sigmaI))
                    *exp(-((y*cos(theta) - x*sin(theta)) * (y*cos(theta) - x*sin(theta)))
                    /(2*sigmaJ*sigmaJ)));
    }

private:
    float sigmaI, sigmaJ, theta;
};


//A 2 dimensional gaussian with spatially varying standard deviations represented
//by polynomials sigmaI(x, y) and sigmaJ(x, y) in its i and j dimensions respectively.  The 
//guassian's i, j coordinate system is rotated by spatialTheta(x, y) (a polynomial) radians 
//with respect to the image's x, y coordinate system.
//The guassian is not normalized.
class spatiallyVaryingGaussian2D: public gaussian2D{
public:
    spatiallyVaryingGaussian2D(polynomial spatialSigmaI_, polynomial spatialSigmaJ_, 
        polynomial spatialTheta_)
            : spatialSigmaI(spatialSigmaI_), spatialSigmaJ(spatialSigmaJ_)
            , spatialTheta(spatialTheta_){}

    Halide::Expr operator()(Halide::Var x, Halide::Var y){
        return (exp(-((x*cos(spatialTheta(x, y)) + y*sin(spatialTheta(x, y))) *
                    (x*cos(spatialTheta(x, y)) + y*sin(spatialTheta(x, y)))) / 
                    (2*spatialSigmaI(x, y)*spatialSigmaI(x, y)))
                    * exp(-((y*cos(spatialTheta(x, y)) - x*sin(spatialTheta(x, y))) * 
                    (y*cos(spatialTheta(x, y)) - x*sin(spatialTheta(x, y)))) / 
                    (2*spatialSigmaJ(x, y)*spatialSigmaJ(x, y))));
    }

private:
    polynomial spatialSigmaI, spatialSigmaJ, spatialTheta;
};

//Represents a 1 dimensional gaussian with standard deviation sigma.
//The guassian is not normalized.
class gaussian1D{
public:
    gaussian1D(float sigma_): sigma(sigma_) {}
    gaussian1D(){} //empty constructor for use by child

    Halide::Expr operator()(Halide::Var i){
        return exp(-(i*i)/(2*sigma*sigma));
    }
private:
    float sigma;
};

//A 1 dimensional gaussian whose standard deviation is represented by the
//spatially varying polynomial spatialSigma(x, y)
class spatiallyVaryingGaussian1D : public gaussian1D{
public:
    spatiallyVaryingGaussian1D(polynomial spatialSigma_): spatialSigma(spatialSigma_) {}

    Halide::Expr operator()(Halide::Var x){
        return exp(-(i*i)/(2*spatialSigma(x, y)*spatialSigma(x, y)));
    }
private:
    polynomial spatialSigma;
};

#endif //__GENERAL_KERNEL_H__
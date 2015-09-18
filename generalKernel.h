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

#define PI_FLOAT 3.14159265359f //
#define NUMBER_OF_RUNS 5 //number of runs when performance testing

//Each output mask pixel is the bitwise OR of each input mask pixel that
//overlaps a nonzero kernel element for the given output pixel.
//If OR_ALL_MASK_PIXELS is true, all mask pixels overlapping a kernel 
//element will be ORed without checking whether the kernel element is zero.
//Gaussians only reach zero at infinity, so this should technically be ok.
//The LSST implementation appears to calculate the kernel with doubles
//and OR when the double is nonzero, but this cutoff seems arbitrary.
//Look into an appropriate cutoff value close to zero.
#define OR_ALL_MASK_PIXELS true


// from lesson_12_using_the_gpu.cpp
// We're going to want to schedule a pipeline in several ways, so we
// define the pipeline in a class so that we can recreate it several
// times with different schedules.

// Define some Vars to use.
Var x, y, i, j, i_v, y_0, yi;



//represents polynomial:
//polynomial(x, y) = a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y 
//                    uvv*x*y*y + vvv*y*y*y
class polynomial{
public:
    polynomial(float a_, float u_, float v_, float uu_, float uv_, float vv_, float uuu_,
        float uuv_, float uvv_, float vvv_): a(a_),  u(u_),  v(v_),  uu(uu_),  
        uv(uv_), vv(vv_),  uuu(uuu_), uuv(uuv_),  uvv(uvv_),  vvv(vvv_) {}

    Halide::Expr operator()(){
        return (a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y +
                    uvv*x*y*y + vvv*y*y*y);
    }

private:
    float a, u, v, uu, uv, vv, uuu, uuv, uvv, vvv;

};

//Abstract base class for 2 dimensional kernels
class kernel2D{
public:
    //return a Halide Expr for the kernel at location (i, j) in the kernel
    virtual Halide::Expr operator()(int i, int j) = 0;

    //return a Halide Expr for the kernel in terms of Vars i and j
    //used for interpolation
    virtual Halide::Expr operator()() = 0;

};



//Represents a 2 dimensional gaussian with standard deviations sigmaI and sigmaJ in its
//i and j dimensions respectively.  The guassian's i, j coordinate system is rotated by theta 
//radians with respect to the image's x, y coordinate system.
//The guassian is normalized from -infinity to infinity, but will not be normalized for 
//a finite kernel size.  Normalization could be removed to save dividing by the normalization
//factor, but is here so that output matches the LSST reference more closely (this came up
//as an issue when perfoming linear interpolation with a spatially varying gaussian)
class gaussian2D: public kernel2D{
public:
    gaussian2D(float sigmaI_, float sigmaJ_, float theta_)
        : sigmaI(sigmaI_), sigmaJ(sigmaJ_), theta(theta_){}

    Halide::Expr operator()(int i, int j){
        return(exp(-((i*cos(theta) + j*sin(theta))*(i*cos(theta) + j*sin(theta)))
                    /(2*sigmaI*sigmaI))
                    *exp(-((j*cos(theta) - i*sin(theta))*(j*cos(theta) - i*sin(theta)))
                    /(2*sigmaJ*sigmaJ)) / (2.0f*PI_FLOAT*sigmaI*sigmaJ));
    }

    Halide::Expr operator()(){
        return(exp(-((i*cos(theta) + j*sin(theta))*(i*cos(theta) + j*sin(theta)))
                    /(2*sigmaI*sigmaI))
                    *exp(-((j*cos(theta) - i*sin(theta))*(j*cos(theta) - i*sin(theta)))
                    /(2*sigmaJ*sigmaJ)) / (2.0f*PI_FLOAT*sigmaI*sigmaJ));
    }

private:
    float sigmaI, sigmaJ, theta;
};


//A 2 dimensional gaussian with spatially varying standard deviations represented
//by polynomials sigmaI(x, y) and sigmaJ(x, y) in its i and j dimensions respectively.  The 
//guassian's i, j coordinate system is rotated by spatialTheta(x, y) (a polynomial) radians 
//with respect to the image's x, y coordinate system.
//The guassian is normalized from -infinity to infinity, but will not be normalized for 
//a finite kernel size.  Normalization could be removed to save dividing by the normalization
//factor, but is here so that output matches the LSST reference more closely (this came up
//as an issue when perfoming linear interpolation with a spatially varying gaussian)
class spatiallyVaryingGaussian2D: public kernel2D{
public:
    spatiallyVaryingGaussian2D(polynomial spatialSigmaI_, polynomial spatialSigmaJ_, 
        polynomial spatialTheta_)
            : spatialSigmaI(spatialSigmaI_), spatialSigmaJ(spatialSigmaJ_)
            , spatialTheta(spatialTheta_){}

/*    Halide::Expr operator()(int i, int j){
        return (exp(-((i*cos(spatialTheta()) + j*sin(spatialTheta())) *
                    (i*cos(spatialTheta()) + j*sin(spatialTheta()))) / 
                    (2*spatialSigmaI()*spatialSigmaI()))
                    * exp(-((j*cos(spatialTheta()) - i*sin(spatialTheta())) * 
                    (j*cos(spatialTheta()) - i*sin(spatialTheta()))) / 
                    (2*spatialSigmaJ()*spatialSigmaJ())))
                    /(2.0f*PI_FLOAT*spatialSigmaI()*spatialSigmaJ());
    }

    Halide::Expr operator()(){
        return (exp(-((i*cos(spatialTheta()) + j*sin(spatialTheta())) *
                    (i*cos(spatialTheta()) + j*sin(spatialTheta()))) / 
                    (2*spatialSigmaI()*spatialSigmaI()))
                    * exp(-((j*cos(spatialTheta()) - i*sin(spatialTheta())) * 
                    (j*cos(spatialTheta()) - i*sin(spatialTheta()))) / 
                    (2*spatialSigmaJ()*spatialSigmaJ())))
                    /(2.0f*PI_FLOAT*spatialSigmaI()*spatialSigmaJ());
    }
*/


    Halide::Expr operator()(int i, int j){
        return (exp(-((i*cos(spatialTheta()) + j*sin(spatialTheta())) *
                    (i*cos(spatialTheta()) + j*sin(spatialTheta()))) / 
                    (2*spatialSigmaI()*spatialSigmaI()))
                    * exp(-((j*cos(spatialTheta()) - i*sin(spatialTheta())) * 
                    (j*cos(spatialTheta()) - i*sin(spatialTheta()))) / 
                    (2*spatialSigmaJ()*spatialSigmaJ())));
    }

    Halide::Expr operator()(){
        return (exp(-((i*cos(spatialTheta()) + j*sin(spatialTheta())) *
                    (i*cos(spatialTheta()) + j*sin(spatialTheta()))) / 
                    (2*spatialSigmaI()*spatialSigmaI()))
                    * exp(-((j*cos(spatialTheta()) - i*sin(spatialTheta())) * 
                    (j*cos(spatialTheta()) - i*sin(spatialTheta()))) / 
                    (2*spatialSigmaJ()*spatialSigmaJ())));
    }

private:
    polynomial spatialSigmaI, spatialSigmaJ, spatialTheta;
};

//Abstract base class for 1 dimensional kernels
class kernel1D{
public:
    virtual Halide::Expr operator()(int i) = 0;
};

//Represents a 1 dimensional gaussian with standard deviation sigma.
//The guassian is not normalized.
class gaussian1D:public kernel1D{
public:
    gaussian1D(float sigma_): sigma(sigma_) {}

    Halide::Expr operator()(int i){
        return exp(-(i*i)/(2*sigma*sigma));
    }
private:
    float sigma;
};

//A 1 dimensional gaussian whose standard deviation is represented by the
//spatially varying polynomial spatialSigma(x, y)
class spatiallyVaryingGaussian1D : public kernel1D{
public:
    spatiallyVaryingGaussian1D(polynomial spatialSigma_): spatialSigma(spatialSigma_) {}

    Halide::Expr operator()(int i){
        return exp(-(i*i)/(2*spatialSigma()*spatialSigma()));
    }
private:
    polynomial spatialSigma;
};



//the total kernel is computed as a (possibly spatially varying) linear combination of each basis
//kernel (if there are multiple basis kernels) and convolved once with each input plane.
//Using a tuple to combine the computation of the three input planes is optional.  Tuples are
//generally faster on the CPU but slower on the GPU.
class generalKernel {
public:
    Image<float> image;
    Image<float> variance;
    Image<uint16_t> mask;

    //output planes calculated using test_performance_cpu()
    Image<float> image_output_cpu;
    Image<float> variance_output_cpu;
    Image<uint16_t> mask_output_cpu;

    //output planes calculated using test_performance_gpu()
    Image<float> image_output_gpu;
    Image<float> variance_output_gpu;
    Image<uint16_t> mask_output_gpu;

    generalKernel(string imageLocation, int kernelSize);

    //Functions to create specific types of kernels

    virtual void createLinearCombinationProgram(
        vector<polynomial> weights, vector<kernel2D *> kernels);

    virtual void createLinearCombinationProgramWithInterpolation(
        vector<polynomial> weights, vector<kernel2D *> kernels, int interpDist);

    virtual void createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(
        vector<polynomial> weights, vector<kernel2D *> kernels, int interpDist);

    //Save .fits images using the LSST stack (does nothing if STANDALONE defined)
    //Before running, load the LSST stack using
    //$ source ./loadLSST.bash
    //$ setup afw -t v10_1   (or appropriate version, also may work with simply $setup afw)
    //save_fits_image("folder1/folder2/imageName") will save image to:
    //folder1/folder2/imageNameKERNELSIZExKERNELSIZE.fits
    //e.g. for the case of a 5x5 kernel:
    //folder1/folder2/imageName5x5.fits
    void save_fits_image(string imageDestination);

    //Kernel schedules
    virtual void schedule_for_cpu() = 0;
    virtual void schedule_for_gpu() = 0;

    //Benchmark kernels and get output image
    virtual void test_performance_cpu() = 0;
    virtual void test_performance_gpu() = 0;

    //check whether the CPU output matches the .fits image
    //referred to by referenceLocation
    //enter 0 for little information, 1 for more details
    void test_correctness(string referenceLocation, int details);

    virtual void debug();

protected:
    //Kernel is size (bounding_box*2 + 1) x (bounding_box*2 + 1)
    int bounding_box;

    //Functions used to bound input planes
    Func image_bounded, variance_bounded, mask_bounded;

    //For interpolation
    Func compressedKernel;

    //Halide expressions used to compute output functions
    Expr total_image_output;
    Expr total_variance_output;
    Expr total_mask_output;
};

class generalKernelWithTuples: public generalKernel{
public:
    generalKernelWithTuples(string imageLocation, int kernelSize)
        : generalKernel(imageLocation, kernelSize){}

    void createLinearCombinationProgram(vector<polynomial> weights, vector<kernel2D *> kernels);
    void createLinearCombinationProgramWithInterpolation(
        vector<polynomial> weights, vector<kernel2D *> kernels, int interpDist);
    void createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(
        vector<polynomial> weights, vector<kernel2D *> kernels, int interpDist);


    void schedule_for_cpu();
    void schedule_for_gpu();

    void schedule_interpolation_for_cpu();


    void test_performance_cpu();
    void test_performance_gpu();

    void debug();

protected:
    //Tuple containing final output of all 3 planes
    Func combined_output;

};

class generalKernelWithoutTuples: public generalKernel{
public:
    generalKernelWithoutTuples(string imageLocation, int kernelSize)
        : generalKernel(imageLocation, kernelSize){}

    void createLinearCombinationProgram(vector<polynomial> weights, vector<kernel2D *> kernels);
    void createLinearCombinationProgramWithInterpolation(
        vector<polynomial> weights, vector<kernel2D *> kernels, int interpDist);
    void createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(
        vector<polynomial> weights, vector<kernel2D *> kernels, int interpDist);


    void schedule_for_cpu();
    void schedule_for_gpu();

    void schedule_interpolation_for_cpu();


    void test_performance_cpu();
    void test_performance_gpu();

    void debug();

protected:
    //Final outputs of the planes without using a tuple
    Func image_output, variance_output, mask_output;
};




#endif //__GENERAL_KERNEL_H__
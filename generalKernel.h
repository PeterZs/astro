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

// from lesson_12_using_the_gpu.cpp
// We're going to want to schedule a pipeline in several ways, so we
// define the pipeline in a class so that we can recreate it several
// times with different schedules.

// Define some Vars to use.
Var x, y, i_v, y_0, yi;



//represents polynomial:
//polynomial(x, y) = a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y 
//                    uvv*x*y*y + vvv*y*y*y
class polynomial{
public:
    polynomial(float a_, float u_, float v_, float uu_, float uv_, float vv_, float uuu_,
        float uuv_, float uvv_, float vvv_): a(a_),  u(u_),  v(v_),  uu(uu_),  
        uv(uv_), vv(vv_),  uuu(uuu_), uuv(uuv_),  uvv(uvv_),  vvv(vvv_) {}

//    Halide::Expr operator()(Halide::Var x, Halide::Var y){
    Halide::Expr operator()(){
        return (a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y +
                    uvv*x*y*y + vvv*y*y*y);
    }


private:
    float a, u, v, uu, uv, vv, uuu, uuv, uvv, vvv;

};

//returns the third degree polynomial given by the input coefficients as a Halide Expr
Halide::Expr getHalidePolynomial(float a, float u, float v, float uu, float uv, float vv,
    float uuu, float uuv, float uvv, float vvv){
    return (a + u*x + v*y + uu*x*x + uv*x*y + vv*y*y + uuu*x*x*x + uuv*x*x*y +uvv*x*y*y +
        vvv*y*y*y);
}

//returns the third degree polynomial with coefficients given by the input vector
//as a Halide Expr
Halide::Expr getHalidePolynomial(vector<float> coefficients){
    if(coefficients.size() != 10)
        cout << "Error getHalidePolynomial, coefficient vector does not have 10 coefficients" << endl;
    return (coefficients[0] + coefficients[1]*x + coefficients[2]*y + coefficients[3]*x*x +
    coefficients[4]*x*y + coefficients[5]*y*y + coefficients[6]*x*x*x + coefficients[7]*x*x*y +
    coefficients[8]*x*y*y + coefficients[9]*y*y*y);
}


//Represents a 2 dimensional gaussian with standard deviations sigmaI and sigmaJ in its
//i and j dimensions respectively.  The guassian's i, j coordinate system is rotated by theta 
//radians with respect to the image's x, y coordinate system.
//The guassian is not normalized.
class gaussian2D{
public:
    gaussian2D(float sigmaI_, float sigmaJ_, float theta_)
        : sigmaI(sigmaI_), sigmaJ(sigmaJ_), theta(theta_){}

    gaussian2D(){} //empty constructor for use by child

//    Halide::Expr operator()(Halide::Var x, Halide::Var y){
    Halide::Expr operator()(int i, int j){
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
//The guassian is not normalized.
class spatiallyVaryingGaussian2D: public gaussian2D{
public:
    spatiallyVaryingGaussian2D(polynomial spatialSigmaI_, polynomial spatialSigmaJ_, 
        polynomial spatialTheta_)
            : spatialSigmaI(spatialSigmaI_), spatialSigmaJ(spatialSigmaJ_)
            , spatialTheta(spatialTheta_){}

//    Halide::Expr operator()(Halide::Var x, Halide::Var y, Halide::Var i, Halide::Var j){
    Halide::Expr operator()(int i, int j){
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

//Represents a 1 dimensional gaussian with standard deviation sigma.
//The guassian is not normalized.
class gaussian1D{
public:
    gaussian1D(float sigma_): sigma(sigma_) {}
    gaussian1D(){} //empty constructor for use by child

//    Halide::Expr operator()(Halide::Var i){
    Halide::Expr operator()(int i){
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

//    Halide::Expr operator()(Halide::Var x, Halide::Var y, Halide::Var i){
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
//    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, kernel1, kernel2,
//    kernel3, kernel4, kernel5, blurImage1, blurImage2, blurImage3, blurImage4, blurImage5,
//    combined_output, image_output, variance_output,
//    mask_output;

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

    //for tuples/no tuples derived classes to call in their constructor 
    generalKernel(){}

    //Functions to create specific types of kernels



    virtual void createLinearCombinationProgram_POOR_NORMALIZATION(
        vector<polynomial> weights, vector<gaussian2D> kernels) = 0;

    virtual void createLinearCombinationProgram_POOR_NORMALIZATION_FAST(
        vector<polynomial> weights, vector<gaussian2D> kernels) = 0;

    virtual void createLinearCombinationProgram(
        vector<polynomial> weights, vector<gaussian2D> kernels) = 0;

    //Save .fits images using the LSST stack
    //Before running, load the LSST stack using
    //$ source ./loadLSST.bash
    //$ setup afw -t v10_1   (or appropriate version, also may work with simply $setup afw)
    void save_fits_image(string imageDestination);

    //Kernel schedules
    virtual void schedule_for_cpu() = 0;
	virtual void schedule_for_gpu() = 0;

    //Benchmark kernels and get output image
    virtual void test_performance_cpu() = 0;
    virtual void test_performance_gpu() = 0;

    //check whether the image planes match
    //implement more if desired later
    void test_correctness(Image<float> reference_output);

    void debug();

protected:
    //Kernel is size (bounding_box*2 + 1) x (bounding_box*2 + 1)
    int bounding_box;

    //Functions used to bound input planes
    Func image_bounded, variance_bounded, mask_bounded;

    //Halide expressions used to compute output functions
    Expr total_image_output;
    Expr total_variance_output;
    Expr total_mask_output;
};

class generalKernelWithTuples: public generalKernel{
public:
    generalKernelWithTuples(string imageLocation, int kernelSize)
        : generalKernel(imageLocation, kernelSize){}

    void createLinearCombinationProgram_VERBOSE(
        vector<polynomial> weights, vector<gaussian2D> kernels);


    void createLinearCombinationProgram_POOR_NORMALIZATION(
        vector<polynomial> weights, vector<gaussian2D> kernels);

    void createLinearCombinationProgram_POOR_NORMALIZATION_FAST(
        vector<polynomial> weights, vector<gaussian2D> kernels);

    void createLinearCombinationProgram(vector<polynomial> weights, vector<gaussian2D> kernels);

    void schedule_for_cpu();
    void schedule_for_gpu();
//
    void test_performance_cpu();
    void test_performance_gpu();

protected:
    //Tuple containing final output of all 3 planes
    Func combined_output;

};

class generalKernelWithoutTuples: public generalKernel{
public:
    generalKernelWithoutTuples(string imageLocation, int kernelSize)
        : generalKernel(imageLocation, kernelSize){}

    void createLinearCombinationProgram(vector<polynomial> weights, vector<gaussian2D> kernels);

    void schedule_for_cpu();
    void schedule_for_gpu();
//
    void test_performance_cpu();
    void test_performance_gpu();

protected:
    //Final outputs of the planes without using a tuple
    Func image_output, variance_output, mask_output;
};




#endif //__GENERAL_KERNEL_H__
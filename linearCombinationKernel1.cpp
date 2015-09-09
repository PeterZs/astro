//for testing without output image 
//compile using: $make blurFast.standalone1
//run using: $./blurFast.standalone1

// On os x:
// g++ linearCombinationKernel1.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o linearCombinationKernel1 -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./linearCombinationKernel1
//
// On linux:
// g++ linearCombinationKernel1.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o linearCombinationKernel1 -std=c++11
//
//LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./linearCombinationKernel1

//this kernel is a spatially varying linear combination of guassians that 
//uses tuples for fast evaluation

#include "linearCombinationKernel1.h"


convolveKernelsSeparatelyThenCombinePipeline::convolveKernelsSeparatelyThenCombinePipeline(Image<float> image_, Image<float> variance_,
    Image<uint16_t> mask_): image(image_), variance(variance_), mask(mask_) {

    //Polynomials that define weights of spatially varying linear combination of 5 kernels
    polynomial1(x, y) = 0.1f + 0.002f*x + 0.003f*y + 0.4f*x*x + 0.5f*x*y
                     + 0.6f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
                     + 0.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial2(x, y) = 1.1f + 1.002f*x + 1.003f*y + 1.4f*x*x + 1.5f*x*y
                     + 1.6f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
                     + 1.00011f*y*y*y;

    //for experimenting with optimizations

    polynomial3(x, y) = 2.1f + 2.002f*x + 2.003f*y + 2.4f*x*x + 2.5f*x*y
                     + 2.6f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
                     + 2.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial4(x, y) = 3.1f + 3.002f*x + 3.003f*y + 3.4f*x*x + 3.5f*x*y
                     + 3.6f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
                     + 3.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial5(x, y) = 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
                     + 4.00011f*y*y*y;

    //5 Kernels that will be weighted by their corresponding polynomials to produce
    //the total kernel
    //Kernel #1
    float sigmaX1 = 2.0f;
    float sigmaY1 = 2.0f;
    float theta1 = 0.0f; //rotation of sigmaX axis
    kernel1(i, j) = (exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
                    *(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1));

    
    //Kernel #2
    float sigmaX2 = 0.5f;
    float sigmaY2 = 4.0f;
    float theta2 = 0.0f; //rotation of sigmaX axis
    kernel2(i, j) = (exp(-((i*cos(theta2) +j*sin(theta2))*(i*cos(theta2) +j*sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (sqrtf(2*M_PI)*sigmaX2))
                    *(exp(-((j*cos(theta2) - i*sin(theta2))*(j*cos(theta2) - i*sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (sqrtf(2*M_PI)*sigmaY2));

    //Kernel #3
    float sigmaX3 = 0.5f;
    float sigmaY3 = 4.0f;
    float theta3 = 3.14159f/4; //rotation of sigmaX axis
    kernel3(i, j) = (exp(-((i*cos(theta3) +j*sin(theta3))*(i*cos(theta3) +j*sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (sqrtf(2*M_PI)*sigmaX3))
                    *(exp(-((j*cos(theta3) - i*sin(theta3))*(j*cos(theta3) - i*sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (sqrtf(2*M_PI)*sigmaY3));
    //Kernel #4
    float sigmaX4 = 0.5f;
    float sigmaY4 = 4.0f;
    float theta4 = 3.14159f/2; //rotation of sigmaX axis
    kernel4(i, j) = (exp(-((i*cos(theta4) +j*sin(theta4))*(i*cos(theta4) +j*sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (sqrtf(2*M_PI)*sigmaX4))
                    *(exp(-((j*cos(theta4) - i*sin(theta4))*(j*cos(theta4) - i*sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (sqrtf(2*M_PI)*sigmaY4));


    //Kernel #5
    float sigmaX5 = 4.0f;
    float sigmaY5 = 4.0f;
    float theta5 = 0.0; //rotation of sigmaX axis
    kernel5(i, j) = (exp(-((i*cos(theta5) +j*sin(theta5))*(i*cos(theta5) +j*sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (sqrtf(2*M_PI)*sigmaX5))
                    *(exp(-((j*cos(theta5) - i*sin(theta5))*(j*cos(theta5) - i*sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (sqrtf(2*M_PI)*sigmaY5));


    //Compute output image plane
    image_bounded = BoundaryConditions::repeat_edge(image);

    //Compute the convolution of each spatially invariant kernel with the image plane
    Expr blur_image_help1 = 0.0f;
    Expr norm1 = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_image_help1 += image_bounded(x + i, y + j) * kernel1(i, j); 
            norm1 += kernel1(i, j);
        }
    }
    blur_image_help1 = blur_image_help1/norm1;
    blurImage1(x, y) = blur_image_help1;

    Expr blur_image_help2 = 0.0f;
    Expr norm2 = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_image_help2 += image_bounded(x + i, y + j) * kernel2(i, j); 
            norm2 += kernel2(i, j);
        }
    }
    blur_image_help2 = blur_image_help2/norm2;
    blurImage2(x, y) = blur_image_help2;

    Expr blur_image_help3 = 0.0f;
    Expr norm3 = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_image_help3 += image_bounded(x + i, y + j) * kernel3(i, j); 
            norm3 += kernel3(i, j);
        }
    }
    blur_image_help3 = blur_image_help3/norm3;
    blurImage3(x, y) = blur_image_help3;

    Expr blur_image_help4 = 0.0f;
    Expr norm4 = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_image_help4 += image_bounded(x + i, y + j) * kernel4(i, j); 
            norm4 += kernel4(i, j);
        }
    }
    blur_image_help4 = blur_image_help4/norm4;
    blurImage4(x, y) = blur_image_help4;

    Expr blur_image_help5 = 0.0f;
    Expr norm5 = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_image_help5 += image_bounded(x + i, y + j) * kernel5(i, j); 
            norm5 += kernel5(i, j);
        }
    }
    blur_image_help5 = blur_image_help5/norm5;
    blurImage5(x, y) = blur_image_help5;
    

    //compute spatially variant linear combination of 5 convolved images
    Expr blur_image_help = (blurImage1(x, y)*polynomial1(x, y) + blurImage2(x, y)*polynomial2(x, y)
                     + blurImage3(x, y)*polynomial3(x, y) + blurImage4(x, y)*polynomial4(x, y) 
                     + blurImage5(x, y)*polynomial5(x, y))/ (polynomial1(x, y) 
                        + polynomial2(x, y) + polynomial3(x, y) + polynomial4(x, y)
                        + polynomial5(x, y));

    //for test speed
//    Expr blur_image_help = (blurImage1(x, y) + blurImage2(x, y)
//                     + blurImage3(x, y) + blurImage4(x, y) 
//                     + blurImage5(x, y))/(polynomial1(x, y) 
//                        + polynomial2(x, y) + polynomial3(x, y) + polynomial4(x, y)
//                        + polynomial5(x, y));

    //Write real variance, mask, if this looks promising
    Expr fakeVar = variance_bounded(x, y) + 2;
    Expr fakeMask = mask_bounded(x, y) + 2;
//    combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, maskOutHelp);
    combined_output(x, y) = Tuple(blur_image_help, fakeVar, fakeMask);
}

void convolveKernelsSeparatelyThenCombinePipeline::debug(){

    //Check out what is happening
    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
    combined_output.compile_to_lowered_stmt("linearCombinationKernel1BlurImage.html", {image}, HTML);

}

void convolveKernelsSeparatelyThenCombinePipeline::schedule_for_cpu() {
    // Split the y coordinate of the consumer into strips of 4 scanlines:
    combined_output.split(y, y_0, yi, 4);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y_0);
    // Vectorize across x by a factor of four.
    combined_output.vectorize(x, 8);

//      blurImage1.compute_root();
//      blurImage2.compute_root();
//      blurImage3.compute_root();
//      blurImage4.compute_root();
//      blurImage5.compute_root();
}

// Now a schedule that uses CUDA or OpenCL.
void convolveKernelsSeparatelyThenCombinePipeline::schedule_for_gpu() {
    // Compute curved in 2D 8x8 tiles using the GPU.
    combined_output.gpu_tile(x, y, 8, 8);

    // JIT-compile the pipeline for the GPU. CUDA or OpenCL are
    // not enabled by default. We have to construct a Target
    // object, enable one of them, and then pass that target
    // object to compile_jit. Otherwise your CPU will very slowly
    // pretend it's a GPU, and use one thread per output pixel.

    // Start with a target suitable for the machine you're running
    // this on.
    Target target = get_host_target();

    // Then enable OpenCL or CUDA.

    // We'll enable OpenCL here, because it tends to give better
    // performance than CUDA, even with NVidia's drivers, because
    // NVidia's open source LLVM backend doesn't seem to do all
    // the same optimizations their proprietary compiler does.
    target.set_feature(Target::OpenCL);

    // Uncomment the next line and comment out the line above to
    // try CUDA instead.
    // target.set_feature(Target::CUDA);

    // If you want to see all of the OpenCL or CUDA API calls done
    // by the pipeline, you can also enable the Debug
    // flag. This is helpful for figuring out which stages are
    // slow, or when CPU -> GPU copies happen. It hurts
    // performance though, so we'll leave it commented out.
    // target.set_feature(Target::Debug);

    combined_output.compile_jit(target);
}

void convolveKernelsSeparatelyThenCombinePipeline::test_performance_cpu() {
    // Benchmark the pipeline.
    image_output(image.width(), image.height());
    variance_output(variance.width(), variance.height());
    mask_output(mask.width(), mask.height());

    Realization r = combined_output.realize(image.width(), image.height());
    image_output = r[0];
    variance_output = r[1];
    mask_output = r[2];

    double average = 0;
    double min;
    double max;
    double imgTime;
    double varTime;
    double maskTime;
    for (int i = 0; i < NUMBER_OF_RUNS; i++) {
        double t1 = current_time();
        r = combined_output.realize(image.width(), image.height());
        double t2 = current_time();
        double t3 = current_time();
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
    average = average/NUMBER_OF_RUNS;
    std::cout << "Average Time: " << average << ", Min = " <<
    min << ", Max = " << max << ", with " << NUMBER_OF_RUNS <<
    " runs" << '\n';
    cout << "For fastest run total time = " << min << ", imgTime = " << imgTime << ", varTime = " << varTime << 
    "maskTime = " << maskTime << endl;

/*
        // If we realize curved into a Halide::Image, that will
        // unfairly penalize GPU performance by including a GPU->CPU
        // copy in every run. Halide::Image objects always exist on
        // the CPU.

        // Halide::Buffer, however, represents a buffer that may
        // exist on either CPU or GPU or both.
        Buffer output(UInt(8), input.width(), input.height(), input.channels());

        // Run the filter once to initialize any GPU runtime state.
        curved.realize(output);

        // Now take the best of 3 runs for timing.
        double best_time;
        for (int i = 0; i < 3; i++) {

            double t1 = current_time();

            // Run the filter 100 times.
            for (int j = 0; j < 100; j++) {
                curved.realize(output);
            }

            // Force any GPU code to finish by copying the buffer back to the CPU.
            output.copy_to_host();

            double t2 = current_time();

            double elapsed = (t2 - t1)/100;
            if (i == 0 || elapsed < best_time) {
                best_time = elapsed;
            }
        }

        printf("%1.4f milliseconds\n", best_time);
*/    
}




convolveOneSpatiallyVaryingKernelPipeline::convolveOneSpatiallyVaryingKernelPipeline(Image<float> image_, Image<float> variance_,
    Image<uint16_t> mask_): image(image_), variance(variance_), mask(mask_) {

    //Polynomials that define weights of spatially variant linear combination of 5 kernels
    polynomial1(x, y) = 0.1f + 0.002f*x + 0.003f*y + 0.4f*x*x + 0.5f*x*y
                     + 0.6f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
                     + 0.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial2(x, y) = 1.1f + 1.002f*x + 1.003f*y + 1.4f*x*x + 1.5f*x*y
                     + 1.6f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
                     + 1.00011f*y*y*y;

    //for experimenting with optimizations

    polynomial3(x, y) = 2.1f + 2.002f*x + 2.003f*y + 2.4f*x*x + 2.5f*x*y
                     + 2.6f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
                     + 2.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial4(x, y) = 3.1f + 3.002f*x + 3.003f*y + 3.4f*x*x + 3.5f*x*y
                     + 3.6f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
                     + 3.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial5(x, y) = 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
                     + 4.00011f*y*y*y;

    //5 Kernels that will be weighted by their corresponding polynomials to produce
    //the total kernel
    //Kernel #1
    float sigmaX1 = 2.0f;
    float sigmaY1 = 2.0f;
    float theta1 = 0.0f; //rotation of sigmaX axis
    kernel1(i, j) = (exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
                    *(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1));

    
    //Kernel #2
    float sigmaX2 = 0.5f;
    float sigmaY2 = 4.0f;
    float theta2 = 0.0f; //rotation of sigmaX axis
    kernel2(i, j) = (exp(-((i*cos(theta2) +j*sin(theta2))*(i*cos(theta2) +j*sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (sqrtf(2*M_PI)*sigmaX2))
                    *(exp(-((j*cos(theta2) - i*sin(theta2))*(j*cos(theta2) - i*sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (sqrtf(2*M_PI)*sigmaY2));

    //Kernel #3
    float sigmaX3 = 0.5f;
    float sigmaY3 = 4.0f;
    float theta3 = 3.14159f/4; //rotation of sigmaX axis
    kernel3(i, j) = (exp(-((i*cos(theta3) +j*sin(theta3))*(i*cos(theta3) +j*sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (sqrtf(2*M_PI)*sigmaX3))
                    *(exp(-((j*cos(theta3) - i*sin(theta3))*(j*cos(theta3) - i*sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (sqrtf(2*M_PI)*sigmaY3));
    //Kernel #4
    float sigmaX4 = 0.5f;
    float sigmaY4 = 4.0f;
    float theta4 = 3.14159f/2; //rotation of sigmaX axis
    kernel4(i, j) = (exp(-((i*cos(theta4) +j*sin(theta4))*(i*cos(theta4) +j*sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (sqrtf(2*M_PI)*sigmaX4))
                    *(exp(-((j*cos(theta4) - i*sin(theta4))*(j*cos(theta4) - i*sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (sqrtf(2*M_PI)*sigmaY4));


    //Kernel #5
    float sigmaX5 = 4.0f;
    float sigmaY5 = 4.0f;
    float theta5 = 0.0; //rotation of sigmaX axis
    kernel5(i, j) = (exp(-((i*cos(theta5) +j*sin(theta5))*(i*cos(theta5) +j*sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (sqrtf(2*M_PI)*sigmaX5))
                    *(exp(-((j*cos(theta5) - i*sin(theta5))*(j*cos(theta5) - i*sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (sqrtf(2*M_PI)*sigmaY5));
    
    
    //Compute output image plane
    image_bounded = BoundaryConditions::repeat_edge(image);    
    Expr blur_image_help = 0.0f;
    Expr norm = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_image_help += image_bounded(x + i, y + j) * (polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j)); 
            norm += (polynomial1(x, y)*kernel1(i, j) + polynomial2(x, y)*kernel2(i, j) + 
                polynomial3(x, y)*kernel3(i, j) + polynomial4(x, y)*kernel4(i, j) + 
                polynomial5(x, y)*kernel5(i, j));
        }
    }
    blur_image_help = blur_image_help/norm;




    //Compute output variance plane
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    Expr blur_variance_help = 0.0f;
//    Expr vNorm2 = 0.0f;
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_variance_help += variance_bounded(x + i, y + j) * (polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j))
                *(polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j)); 
//            vNorm2 += (polynomial1(x, y)*kernel1(i, j) + polynomial2(x, y)*kernel2(i, j) + 
//                polynomial3(x, y)*kernel3(i, j) + polynomial4(x, y)*kernel4(i, j) + 
//                polynomial5(x, y)*kernel5(i, j));
        }
    }
//    blur_variance_help = blur_variance_help/(norm(x,y)*norm(x,y));
    blur_variance_help = blur_variance_help/(norm*norm);
        


    //Compute output mask plane
    mask_bounded = BoundaryConditions::repeat_edge(mask);    
    Expr blur_mask_help = cast<uint16_t>(0);  //make sure blur_mask_help has type uint16
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_mask_help = select((polynomial1(x, y)*kernel1(i, j) + polynomial2(x, y)*kernel2(i, j) + 
                polynomial3(x, y)*kernel3(i, j) + polynomial4(x, y)*kernel4(i, j) + 
                polynomial5(x, y)*kernel5(i, j)) == 0.0f, blur_mask_help, blur_mask_help | mask_bounded(x + i, y + j));
//            blur_mask_help = blur_mask_help | mask_bounded(x + i, y + j);    
        }
    }

    //Evaluate image, mask, and variance planes concurrently using a tuple
    combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, blur_mask_help);
}

void convolveOneSpatiallyVaryingKernelPipeline::schedule_for_cpu() {
    // Split the y coordinate of the consumer into strips of 4 scanlines:
    combined_output.split(y, y_0, yi, 4);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y_0);
    // Vectorize across x by a factor of four.
    combined_output.vectorize(x, 8);
}

// Now a schedule that uses CUDA or OpenCL.
void convolveOneSpatiallyVaryingKernelPipeline::schedule_for_gpu() {
    // Compute curved in 2D 8x8 tiles using the GPU.
    combined_output.gpu_tile(x, y, 8, 8);

    // JIT-compile the pipeline for the GPU. CUDA or OpenCL are
    // not enabled by default. We have to construct a Target
    // object, enable one of them, and then pass that target
    // object to compile_jit. Otherwise your CPU will very slowly
    // pretend it's a GPU, and use one thread per output pixel.

    // Start with a target suitable for the machine you're running
    // this on.
    Target target = get_host_target();

    // Then enable OpenCL or CUDA.

    // We'll enable OpenCL here, because it tends to give better
    // performance than CUDA, even with NVidia's drivers, because
    // NVidia's open source LLVM backend doesn't seem to do all
    // the same optimizations their proprietary compiler does.
    target.set_feature(Target::OpenCL);

    // Uncomment the next line and comment out the line above to
    // try CUDA instead.
    // target.set_feature(Target::CUDA);

    // If you want to see all of the OpenCL or CUDA API calls done
    // by the pipeline, you can also enable the Debug
    // flag. This is helpful for figuring out which stages are
    // slow, or when CPU -> GPU copies happen. It hurts
    // performance though, so we'll leave it commented out.
    // target.set_feature(Target::Debug);

    combined_output.compile_jit(target);
}

void convolveOneSpatiallyVaryingKernelPipeline::test_performance_cpu() {
    // Test the performance of the pipeline.

 
    // Benchmark the pipeline.
    image_output(image.width(), image.height());
    variance_output(variance.width(), variance.height());
    mask_output(mask.width(), mask.height());

    Realization r = combined_output.realize(image.width(), image.height());
    image_output = r[0];
    variance_output = r[1];
    mask_output = r[2];

    double average = 0;
    double min;
    double max;
    double imgTime;
    double varTime;
    double maskTime;
    for (int i = 0; i < NUMBER_OF_RUNS; i++) {
        double t1 = current_time();
        r = combined_output.realize(image.width(), image.height());
        double t2 = current_time();
        double t3 = current_time();
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
    average = average/NUMBER_OF_RUNS;
    std::cout << "Average Time: " << average << ", Min = " <<
    min << ", Max = " << max << ", with " << NUMBER_OF_RUNS <<
    " runs" << '\n';
    cout << "For fastest run total time = " << min << ", imgTime = " << imgTime << ", varTime = " << varTime << 
    "maskTime = " << maskTime << endl;   
}

void convolveOneSpatiallyVaryingKernelPipeline::test_performance_gpu() {
    // Test the performance of the pipeline.
    // If we realize curved into a Halide::Image, that will
    // unfairly penalize GPU performance by including a GPU->CPU
    // copy in every run. Halide::Image objects always exist on
    // the CPU.

    // Halide::Buffer, however, represents a buffer that may
    // exist on either CPU or GPU or both.
//    Buffer output(UInt(8), input.width(), input.height(), input.channels());

    // Run the filter once to initialize any GPU runtime state.
    Realization r = combined_output.realize(image.width(), image.height());
    image_output = r[0];
    variance_output = r[1];
    mask_output = r[2];

    // Now take the best of 3 runs for timing.
    double best_time;
    for (int i = 0; i < 3; i++) {

        double t1 = current_time();

        // Run the filter 100 times.
        for (int j = 0; j < 100; j++) {
            combined_output.realize(image.width(), image.height());
        }

        // Force any GPU code to finish by copying the buffer back to the CPU.
//        output.copy_to_host();
        //Does this do the equivalent?
        image_output = r[0];
        variance_output = r[1];
        mask_output = r[2];

        double t2 = current_time();

        double elapsed = (t2 - t1)/100;
        if (i == 0 || elapsed < best_time) {
            best_time = elapsed;
        }
    }

    printf("%1.4f milliseconds\n", best_time);
}

void convolveOneSpatiallyVaryingKernelPipeline::debug(){
    //Check out what is happening
    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
    combined_output.compile_to_lowered_stmt("linearCombinationKernel1BlurImage.html", {image}, HTML);
}



int main(int argc, char *argv[]) {

#ifndef STANDALONE
    auto im = afwImage::MaskedImage<float>("./images/calexp-004207-g3-0123.fits");
    int width = im.getWidth(), height = im.getHeight();
#else
    int width = 2048, height = 1489;
    printf("[no load]");
#endif

    printf("Loaded: %d x %d\n", width, height);
    Image<float> image(width, height);
    Image<float> variance(width, height);
    Image<uint16_t> mask(width, height);

#ifndef STANDALONE
    //Read image in
    for (int y = 0; y < height; y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = im.x_at(0, y);
        for (int x = 0; x < width; x++){
            image(x, y) = (*inPtr).image();
            variance(x, y) = (*inPtr).variance();
            mask(x, y) = (*inPtr).mask();
            inPtr++;
        }
    }
#endif
    convolveOneSpatiallyVaryingKernelPipeline p1(image, variance, mask);

    p1.schedule_for_cpu();
    p1.test_performance_cpu();

//    p1.schedule_for_gpu();
//    p1.test_performance_gpu();

//    p1.debug(); 



#ifndef STANDALONE
    //write image out
    auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
    for (int y = 0; y < imOut.getHeight(); y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);

        for (int x = 0; x < imOut.getWidth(); x++){
            afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
            curPixel(p1.image_output(x, y), p1.mask_output(x, y), p1.variance_output(x, y));
            (*inPtr) = curPixel;
            inPtr++;

        }
    }

    imOut.writeFits("./halideLinearCombination1.fits");
#endif

}




//Other polynomials
/*    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial3 ("polynomial3");
    polynomial3(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial4 ("polynomial4");
    polynomial4(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial5 ("polynomial5");
    polynomial5(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;
*/
    //Testing different polynomials
/*    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 1.1f + 2.001f*x + 0.301f*y + 0.00401f*x*x + 0.034501f*x*y
                     + 0.0023451f*y*y +  0.0234534001f*x*x*x + 0.0234500001f*x*x*y + 5.0300001f*x*y*y
                     + 0.000123412000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 05.1f + 0.0401f*x + 0.4031f*y + 03.06401f*x*x + 04.000001f*x*y
                     + 50.000234001f*y*y +  0.002345340001f*x*x*x + 0.054300001f*x*x*y + 0.03400000001f*x*y*y
                     + 0.0634543001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial3 ("polynomial3");
    polynomial3(x, y) = 3.1f + 0.005431f*x + 0.2345001f*y + 0.2345000001f*x*x + 0.532000001f*x*y
                     + 0.003451f*y*y +  0.0005340001f*x*x*x + 0.023450001f*x*x*y + 235.000000001f*x*y*y
                     + 345.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial4 ("polynomial4");
    polynomial4(x, y) = 30.1f + 0.345001f*x + 0.543001f*y + 0.4567f*x*x + 0.2345000001f*x*y
                     + 0.003453401f*y*y +  0.000657860001f*x*x*x + 0.5342000000001f*x*x*y + 0.2345000000001f*x*y*y
                     + 0.5234000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial5 ("polynomial5");
    polynomial5(x, y) = 6.1f + 34.001f*x + 0.543001f*y + 0.34000001f*x*x + 0.534000001f*x*y
                     + 0.345601f*y*y +  0345.000000001f*x*x*x + 0.053400000001f*x*x*y + 0.0003245000001f*x*y*y
                     + 0.0006345001f*y*y*y;

*/


/*    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 1.1f + 1.001f*x + 1.001f*y + 1.000001f*x*x + 1.000001f*x*y
                     + 1.000001f*y*y +  1.000000001f*x*x*x + 1.000000001f*x*x*y + 1.000000001f*x*y*y
                     + 1.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial3 ("polynomial3");
    polynomial3(x, y) = 2.1f + 2.001f*x + 2.001f*y + 2.000001f*x*x + 2.000001f*x*y
                     + 2.000001f*y*y +  2.000000001f*x*x*x + 2.000000001f*x*x*y + 2.000000001f*x*y*y
                     + 2.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial4 ("polynomial4");
    polynomial4(x, y) = 3.1f + 3.001f*x + 3.001f*y + 3.000001f*x*x + 3.000001f*x*y
                     + 3.000001f*y*y +  3.000000001f*x*x*x + 3.000000001f*x*x*y + 3.000000001f*x*y*y
                     + 3.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial5 ("polynomial5");
    polynomial5(x, y) = 4.1f + 4.001f*x + 4.001f*y + 4.000001f*x*x + 4.000001f*x*y
                     + 4.000001f*y*y +  4.000000001f*x*x*x + 4.000000001f*x*x*y + 4.000000001f*x*y*y
                     + 4.000000001f*y*y*y;
*/


    //compute output image and variance
/*    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.001f*x + 0.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 0.1f + 1.001f*x + 1.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial3 ("polynomial3");
    polynomial3(x, y) = 0.1f + 2.001f*x + 2.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial4 ("polynomial4");
    polynomial4(x, y) = 0.1f + 3.001f*x + 3.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;

    //for experimenting with optimizations
    Func polynomial5 ("polynomial5");
    polynomial5(x, y) = 0.1f + 4.001f*x + 4.001f*y + 0.000001f*x*x + 0.000001f*x*y
                     + 0.000001f*y*y +  0.000000001f*x*x*x + 0.000000001f*x*x*y + 0.000000001f*x*y*y
                     + 0.000000001f*y*y*y;


    Func polynomial1 ("polynomial1");
    polynomial1(x, y) = 0.1f + 0.002f*x + 0.003f*y + 0.000004f*x*x + 0.000005f*x*y
                     + 0.000006f*y*y +  0.000000007f*x*x*x + 0.000000008f*x*x*y + 0.000000009f*x*y*y
                     + 0.0000000011f*y*y*y;

    //for experimenting with optimizations
    Func polynomial2 ("polynomial2");
    polynomial2(x, y) = 1.1f + 1.002f*x + 1.003f*y + 1.000004f*x*x + 1.000005f*x*y
                     + 1.000006f*y*y +  1.000000007f*x*x*x + 1.000000008f*x*x*y + 1.000000009f*x*y*y
                     + 1.0000000011f*y*y*y;

    //for experimenting with optimizations
    Func polynomial3 ("polynomial3");
    polynomial3(x, y) = 2.1f + 2.002f*x + 2.003f*y + 2.000004f*x*x + 2.000005f*x*y
                     + 2.000006f*y*y +  2.000000007f*x*x*x + 2.000000008f*x*x*y + 2.000000009f*x*y*y
                     + 2.0000000011f*y*y*y;

    //for experimenting with optimizations
    Func polynomial4 ("polynomial4");
    polynomial4(x, y) = 3.1f + 3.002f*x + 3.003f*y + 3.000004f*x*x + 3.000005f*x*y
                     + 3.000006f*y*y +  3.000000007f*x*x*x + 3.000000008f*x*x*y + 3.000000009f*x*y*y
                     + 3.0000000011f*y*y*y;

    //for experimenting with optimizations
    Func polynomial5 ("polynomial5");
    polynomial5(x, y) = 4.1f + 4.002f*x + 4.003f*y + 4.000004f*x*x + 4.000005f*x*y
                     + 4.000006f*y*y +  4.000000007f*x*x*x + 4.000000008f*x*x*y + 4.000000009f*x*y*y
                     + 4.0000000011f*y*y*y;
*/

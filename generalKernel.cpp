// On os x:
// g++ generalKernel.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o generalKernel -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./generalKernel
//
// On linux:
// g++ generalKernel.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o generalKernel -std=c++11
//
//LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./generalKernel


//fix this
//for testing without output image 
//compile using: $make blurFast.standalone1
//run using: $./blurFast.standalone1

//this kernel is a spatially varying linear combination of guassians that 
//uses tuples for fast evaluation

#include "generalKernel.h"


convolveKernelsSeparatelyThenCombinePipeline::convolveKernelsSeparatelyThenCombinePipeline(
    Image<float> image_, Image<float> variance_, Image<uint16_t> mask_, bool useTuples_):
    image(image_), variance(variance_), mask(mask_), useTuples(useTuples_) {

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
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    mask_bounded = BoundaryConditions::repeat_edge(mask);

    Expr fakeVar = variance_bounded(x, y) + 2;
    Expr fakeMask = mask_bounded(x, y) + 2;
//    combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, maskOutHelp);

    //Evaluate image, mask, and variance planes concurrently using a tuple
    if(useTuples){
//        combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, blur_mask_help);
        combined_output(x, y) = Tuple(blur_image_help, fakeVar, fakeMask);
    }
    else{
        imageOut(x, y) = blur_image_help;
        varianceOut(x, y) = fakeVar;
        maskOut(x, y) = fakeMask;
//        varianceOut(x, y) = blur_variance_help;
//        maskOut(x, y) = blur_mask_help;
    }
}

void convolveKernelsSeparatelyThenCombinePipeline::debug(){

    //Check out what is happening
    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
    combined_output.compile_to_lowered_stmt("generalKernelBlurImage.html", {image}, HTML);

}

void convolveKernelsSeparatelyThenCombinePipeline::schedule_for_cpu() {
    // Split the y coordinate of the consumer into strips of 4 scanlines:
    combined_output.split(y, y_0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y_0);
    // Vectorize across x by a factor of four.
    combined_output.vectorize(x, 16);

//      blurImage1.compute_root();
//      blurImage2.compute_root();
//      blurImage3.compute_root();
//      blurImage4.compute_root();
//      blurImage5.compute_root();
}

// Now a schedule that uses CUDA or OpenCL.
void convolveKernelsSeparatelyThenCombinePipeline::schedule_for_gpu() {
   // Compute curved in 2D 8x8 tiles using the GPU.

    if(useTuples){
        combined_output.gpu_tile(x, y, 16, 16);
    }

    else{
        blurImage1.compute_root();
        blurImage2.compute_root();
        blurImage3.compute_root();
        blurImage4.compute_root();
        blurImage5.compute_root();

        imageOut.gpu_tile(x, y, 16, 16);
        varianceOut.gpu_tile(x, y, 16, 16);
        maskOut.gpu_tile(x, y, 16, 16);
    }

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

    if(useTuples){
        combined_output.compile_jit(target);
    }
    else{
        imageOut.compile_jit(target);
        varianceOut.compile_jit(target);
        maskOut.compile_jit(target);
    }
}

void convolveKernelsSeparatelyThenCombinePipeline::test_performance_cpu() {
    // Benchmark the pipeline.
    Image<float> image_output(image.width(), image.height());
    Image<float> variance_output(variance.width(), variance.height());
    Image<uint16_t> mask_output(mask.width(), mask.height());

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


void convolveKernelsSeparatelyThenCombinePipeline::test_performance_gpu() {
    // Test the performance of the pipeline.
    // If we realize curved into a Halide::Image, that will
    // unfairly penalize GPU performance by including a GPU->CPU
    // copy in every run. Halide::Image objects always exist on
    // the CPU.

    // Halide::Buffer, however, represents a buffer that may
    // exist on either CPU or GPU or both.
//    Buffer output(UInt(8), input.width(), input.height(), input.channels());

    // Run the filter once to initialize any GPU runtime state.
    Image<float> image_output(image.width(), image.height());
    Image<float> variance_output(variance.width(), variance.height());
    Image<uint16_t> mask_output(mask.width(), mask.height());

    //need to define r outside of if statement, won't compile without assignment
    Func fake;
    fake(x,y) = 0;
    Realization r = fake.realize(1, 1);

    if(useTuples){
        r = combined_output.realize(image.width(), image.height());
        image_output = r[0];
        variance_output = r[1];
        mask_output = r[2];
    }

    else{
        image_gpu_output = Buffer(Float(32), image.width(), image.height());
        variance_gpu_output = Buffer(Float(32), image.width(), image.height());
        mask_gpu_output = Buffer(UInt(16), image.width(), image.height());        
    }

    // Now take the best of 3 runs for timing.
    double best_time;
    double t1;
    double t2;
    double elapsed;
    for (int i = 0; i < 3; i++) {

        if(useTuples){
            t1 = current_time();
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
            t2 = current_time();
        }
        else{
            t1 = current_time();
            // Run the filter 100 times.
            for (int j = 0; j < 100; j++) {
                imageOut.realize(image_gpu_output);
                varianceOut.realize(variance_gpu_output);
                maskOut.realize(mask_gpu_output);
            }
            // Force any GPU code to finish by copying the buffer back to the CPU.
            image_gpu_output.copy_to_host();
            variance_gpu_output.copy_to_host();
            mask_gpu_output.copy_to_host();

            t2 = current_time();

        }
        elapsed = (t2 - t1)/100;
        if (i == 0 || elapsed < best_time) {
            best_time = elapsed;
        }
    }

    printf("%1.4f milliseconds\n", best_time);
}




convolveOncePipeline::convolveOncePipeline(
    Image<float> image_, Image<float> variance_, Image<uint16_t> mask_, bool useTuples_):
    image(image_), variance(variance_), mask(mask_), useTuples(useTuples_) {
    //original LSST example
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
    //Polynomials that define weights of spatially varying linear combination of 5 kernels
    polynomial p1(0.1, 0.002, 0.003, 0.4, 0.5, 0.6,  0.0007, 0.0008, 0.0009, 0.00011);
    polynomial p2(1.1, 1.002, 1.003, 1.4, 1.5, 1.6,  1.0007, 1.0008, 1.0009, 1.00011);
    polynomial p3(2.1, 2.002, 2.003, 2.4, 2.5, 2.6,  2.0007, 2.0008, 2.0009, 2.00011);
    polynomial p4(3.1, 3.002, 3.003, 3.4, 3.5, 3.6,  3.0007, 3.0008, 3.0009, 3.00011);
    polynomial p5(4.1, 4.002, 4.003, 4.4, 4.5, 4.6,  4.0007, 4.0008, 4.0009, 4.00011);
    polynomial1(x, y) = p1(x, y);
    polynomial2(x, y) = p2(x, y);
    polynomial3(x, y) = p3(x, y);
    polynomial4(x, y) = p4(x, y);
    polynomial5(x, y) = p5(x, y);

    //5 Kernels that will be weighted by their corresponding polynomials to produce
    //the total kernel
    //Kernel #1
    gaussian2D k1(2, 2, 0);
    kernel1(i, j) = k1(i, j);

    //Kernel #2
    gaussian2D k2(.5, 4, 0);
    kernel2(i, j) = k2(i, j);

    //Kernel #3
    gaussian2D k3(.5, 4, M_PI/4);
    kernel3(i, j) = k3(i, j);

    //Kernel #4
    gaussian2D k4(.5, 4, M_PI/2);
    kernel4(i, j) = k4(i, j);

    //Kernel #5
    gaussian2D k5(4, 4, 0);
    kernel5(i, j) = k5(i, j);

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
    for(int i = -BOUNDING_BOX; i <= BOUNDING_BOX; i++){
        for(int j = -BOUNDING_BOX; j <= BOUNDING_BOX; j++){
            blur_variance_help += variance_bounded(x + i, y + j) * (polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j))
                *(polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j)); 
        }
    }
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
    if(useTuples){
        combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, blur_mask_help);
    }
    else{
        imageOut(x, y) = blur_image_help;
        varianceOut(x, y) = blur_variance_help;
        maskOut(x, y) = blur_mask_help;
    }
}

void convolveOncePipeline::schedule_for_cpu() {
    if(useTuples){
        // Split the y coordinate of the consumer into strips of 4 scanlines:
        combined_output.split(y, y_0, yi, 32);
        // Compute the strips using a thread pool and a task queue.
        combined_output.parallel(y_0);
        // Vectorize across x by a factor of four.
        combined_output.vectorize(x, 16);
    }
    else{
        imageOut.split(y, y_0, yi, 32);
        imageOut.parallel(y_0);
        imageOut.vectorize(x, 16);

        varianceOut.split(y, y_0, yi, 32);
        varianceOut.parallel(y_0);
        varianceOut.vectorize(x, 16);

        maskOut.split(y, y_0, yi, 32);
        maskOut.parallel(y_0);
        maskOut.vectorize(x, 16);
    }   
}

// Now a schedule that uses CUDA or OpenCL.
void convolveOncePipeline::schedule_for_gpu() {
    // Compute curved in 2D 8x8 tiles using the GPU.

    if(useTuples){
        combined_output.gpu_tile(x, y, 16, 16);
    }

    else{
        imageOut.gpu_tile(x, y, 16, 16);
        varianceOut.gpu_tile(x, y, 16, 16);
        maskOut.gpu_tile(x, y, 16, 16);
    }

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

    if(useTuples){
        combined_output.compile_jit(target);
    }
    else{
        imageOut.compile_jit(target);
        varianceOut.compile_jit(target);
        maskOut.compile_jit(target);
    }
}

void convolveOncePipeline::test_performance_cpu() {
    // Test the performance of the pipeline.

 
    // Benchmark the pipeline.
    Image<float> image_output_(image.width(), image.height());
    image_output = image_output_;

//    Image<float> image_output(image.width(), image.height());
    Image<float> variance_output(variance.width(), variance.height());
    Image<uint16_t> mask_output(mask.width(), mask.height());


    //need to define r outside of if statement, won't compile without assignment
    Func fake;
    fake(x,y) = 0;
    Realization r = fake.realize(1, 1);

    if(useTuples){
        r = combined_output.realize(image.width(), image.height());
        image_output = r[0];
        variance_output = r[1];
        mask_output = r[2];
    }
    else{
        Image<float> fake_output(image.width(), image.height());
        imageOut.realize(fake_output);

        image_output(image.width(), image.height());
        imageOut.realize(image_output);
        varianceOut.realize(variance_output);
        maskOut.realize(mask_output);
    }

    double average = 0;
    double min;
    double max;
    double imgTime;
    double varTime;
    double maskTime;

    double t1;
    double t2;
    double t3;
    double t4;
    double curTime;
    for (int i = 0; i < NUMBER_OF_RUNS; i++) {
        if(useTuples){
            t1 = current_time();
            r = combined_output.realize(image.width(), image.height());
            t2 = current_time();
            t3 = current_time();
            t4 = current_time();
            curTime = (t4-t1);
        }
        else{
            t1 = current_time();
            imageOut.realize(image_output);
            t2 = current_time();
            varianceOut.realize(variance_output);
            t3 = current_time();
            maskOut.realize(mask_output);
            t4 = current_time();
            curTime = (t4-t1);
        }
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

void convolveOncePipeline::test_performance_gpu() {
    //initialize output images with correct size
    Image<float> image_output_(image.width(), image.height());
    Image<float> variance_output_(variance.width(), variance.height());
    Image<uint16_t> mask_output_(mask.width(), mask.height());

    image_output = image_output_;
    variance_output = variance_output_;
    mask_output = mask_output_;

    // Test the performance of the pipeline.
    // If we realize curved into a Halide::Image, that will
    // unfairly penalize GPU performance by including a GPU->CPU
    // copy in every run. Halide::Image objects always exist on
    // the CPU.

    // Halide::Buffer, however, represents a buffer that may
    // exist on either CPU or GPU or both.
    Buffer image_gpu_output = Buffer(Float(32), image.width(), image.height());
    Buffer variance_gpu_output = Buffer(Float(32), image.width(), image.height());
    Buffer mask_gpu_output = Buffer(UInt(16), image.width(), image.height());    

    //need to define r outside of if statement, won't compile without assignment
    Func fake;
    fake(x,y) = 0;
    Realization r = fake.realize(1, 1);

    // Run the filter once to initialize any GPU runtime state.
    // Also save the output image, variance, and mask planes
    if(useTuples){
        r = combined_output.realize(image.width(), image.height());
        image_output = r[0];
        variance_output = r[1];
        mask_output = r[2];
    }

    else{
        imageOut.realize(image_gpu_output);
        varianceOut.realize(variance_gpu_output);
        maskOut.realize(mask_gpu_output);
        image_output = image_gpu_output;
        variance_output = variance_gpu_output;
        mask_output = mask_gpu_output;
    }

    // Now take the best of 3 runs for timing.
    double best_time;
    double t1;
    double t2;
    double elapsed;
    for (int i = 0; i < 3; i++) {

        if(useTuples){
            t1 = current_time();
            // Run the filter 100 times.
            for (int j = 0; j < 100; j++) {
                combined_output.realize(image.width(), image.height());
            }
            // Force any GPU code to finish by copying the buffer back to the CPU.
            //does the extra step of copying to a buffer account for noticable time?
            image_gpu_output = r[0];
            variance_gpu_output = r[1];
            mask_gpu_output = r[2];
            image_gpu_output.copy_to_host();
            variance_gpu_output.copy_to_host();
            mask_gpu_output.copy_to_host();
            t2 = current_time();
        }
        else{
            t1 = current_time();
            // Run the filter 100 times.
            for (int j = 0; j < 100; j++) {
                imageOut.realize(image_gpu_output);
                varianceOut.realize(variance_gpu_output);
                maskOut.realize(mask_gpu_output);
            }
            // Force any GPU code to finish by copying the buffer back to the CPU.
            image_gpu_output.copy_to_host();
            variance_gpu_output.copy_to_host();
            mask_gpu_output.copy_to_host();

            t2 = current_time();

        }
        elapsed = (t2 - t1)/100;
        if (i == 0 || elapsed < best_time) {
            best_time = elapsed;
        }
    }

    printf("%1.4f milliseconds\n", best_time);
}

void convolveOncePipeline::debug(){
    //Check out what is happening
    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
    combined_output.compile_to_lowered_stmt("generalKernelBlurImage.html", {image}, HTML);
}

void convolveOncePipeline::test_correctness(Image<float> reference_output) {
    Image<float> image_output(image.width(), image.height());
    imageOut.realize(image_output);


    // Check against the reference output.
    for (int y = 0; y < image.height(); y++) {
        for (int x = 0; x < image.width(); x++) {
            if (image_output(x, y) != reference_output(x, y)) {
                printf("Mismatch between output (%f) and "
                       "reference output (%f) at %d, %d\n",
                       image_output(x, y),
                       reference_output(x, y),
                       x, y);
                exit(0);
            }
        }
    }
    cout << "done checking correctness" << endl;

}

int main(int argc, char *argv[]) {

    cout << "Kernel size = " << (BOUNDING_BOX*2 + 1) << " x " << (BOUNDING_BOX*2 + 1) <<endl;

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
//    convolveKernelsSeparatelyThenCombinePipeline p0(image, variance, mask, false);
//    cout << "convolveKernelsSeparatelyThenCombinePipeline On GPU without tuples: " << endl;
//    p0.schedule_for_gpu();
//    p0.test_performance_gpu();

/*    convolveOncePipeline p1(image, variance, mask, true);
    cout << "On CPU with tuples: " << endl;
    p1.schedule_for_cpu();
    p1.test_performance_cpu();

    convolveOncePipeline p2(image, variance, mask, true);
    cout << "On GPU with tuples: " << endl;
    p2.schedule_for_gpu();
    p2.test_performance_gpu();
*/
    convolveOncePipeline p3(image, variance, mask, false);
    cout << "On CPU without tuples: " << endl;
    p3.schedule_for_cpu();
    p3.test_performance_cpu();

//    convolveOncePipeline p4(image, variance, mask, false);
//    cout << "On GPU without tuples: " << endl;
//    p4.schedule_for_gpu();
//    p4.test_performance_gpu();

    //Check GPU:
    //Allocate an image that will store the correct image plane output 
    //calculated on CPU
    Image<float> reference_output_image(image.width(), image.height());
    p3.imageOut.realize(reference_output_image);
    //check it matches GPU calculation
//    p4.test_correctness(reference_output_image);


    //calculate variance/mask planes for writing out image
    Image<float> reference_output_variance(image.width(), image.height());
    Image<uint16_t> reference_output_mask(image.width(), image.height());
    p3.varianceOut.realize(reference_output_variance);
    p3.maskOut.realize(reference_output_mask);

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
            curPixel(reference_output_image(x, y), reference_output_mask(x, y),
            reference_output_variance(x, y));
            (*inPtr) = curPixel;
            inPtr++;

        }
    }
    std::string B_BOX_STRING = std::to_string(BOUNDING_BOX*2 + 1);
    imOut.writeFits("./images/linearCombination/halideLinearCombination" + B_BOX_STRING +
    "x" + B_BOX_STRING + ".fits");
#endif

}




//Other polynomials

/*
    //original LSST example
    Func polynomial1 ("polynomial1");
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

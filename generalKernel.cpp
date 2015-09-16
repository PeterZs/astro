// On os x:
// g++ generalKernel.cpp -g -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o generalKernel -std=c++11
//
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./generalKernel
//
// On linux:
//g++ generalKernel.cpp -g -I ./include -I Linux64/pex_policy/10.1+1/include/ -I Linux64/daf_persistence/10.1+1/include/ -I Linux64/utils/10.1+1/include/ -I Linux64/daf_base/10.1+2/include/ -I Linux64/base/10.1+1/include/ -I Linux64/ndarray/10.1+2/include/ -I Linux64/pex_exceptions/10.1+1/include/ -I Linux64/eigen/3.2.0/include/ -I Linux64/afw/10.1+1/include -L ./bin -L Linux64/afw/10.1+1/lib -L Linux64/daf_base/10.1+2/lib/ -L Linux64/daf_persistence/10.1+1/lib/ -L Linux64/boost/1.55.0.1.lsst2+3/lib/ -L Linux64/wcslib/4.14+7/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -lpthread -ldl -o generalKernel -std=c++11
//
//LD_LIBRARY_PATH=./bin:Linux64/afw/10.1+1/lib/:Linux64/daf_persistence/10.1+1/lib/:Linux64/daf_base/10.1+2/lib/:Linux64/boost/1.55.0.1.lsst2+3/lib/:Linux64/xpa/2.1.15.lsst2/lib/:Linux64/pex_policy/10.1+1/lib/:Linux64/pex_logging/10.1+1/lib/:Linux64/utils/10.1+1/lib/:Linux64/pex_exceptions/10.1+1/lib/:Linux64/base/10.1+1/lib/:Linux64/wcslib/4.14+7/lib/:Linux64/cfitsio/3360.lsst1/lib/:Linux64/gsl/1.16.lsst1/lib/:Linux64/minuit2/5.28.00/lib:Linux64/mysql/5.1.65.lsst2/lib/ ./generalKernel


//fix this
//for testing without output image 
//compile using: $make blurFast.standalone1
//run using: $./blurFast.standalone1

//this kernel is a spatially varying linear combination of guassians that 
//uses tuples for fast evaluation

#include "generalKernel.h"

generalKernel::generalKernel(string imageLocation, int kernelSize) {

    bounding_box = (kernelSize - 1)/2;

    //load input image, variance, and mask planes
    #ifndef STANDALONE
        auto im = afwImage::MaskedImage<float>(imageLocation);
        int width = im.getWidth(), height = im.getHeight();
    #else
        int width = 2048, height = 1489;
        printf("[no load]");
    #endif

        //initialize image, variance and mask planes to correct size
        printf("Loaded: %d x %d\n", width, height);
        Image<float> image_(width, height);
        Image<float> variance_(width, height);
        Image<uint16_t> mask_(width, height);

        image = image_;
        variance = variance_;
        mask = mask_;

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

    image_bounded = BoundaryConditions::repeat_edge(image);    
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    mask_bounded = BoundaryConditions::repeat_edge(mask);    
}

//save fits image (does nothing if STANDALONE defined)
//save_fits_image("folder1/folder2/imageName") will save image to:
//folder1/folder2/imageNameKERNELSIZExKERNELSIZE.fits
//e.g. for the case of a 5x5 kernel:
//folder1/folder2/imageName5x5.fits
void generalKernel::save_fits_image(string imageDestination){
    #ifndef STANDALONE
        //write image out
        auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel,
                        lsst::afw::image::VariancePixel>(image.width(), image.height());

        for (int y = 0; y < imOut.getHeight(); y++) {
            afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, 
            lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);

            for (int x = 0; x < imOut.getWidth(); x++){
                afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> 
                curPixel(image_output_cpu(x, y), mask_output_cpu(x, y), variance_output_cpu(x, y));
                (*inPtr) = curPixel;
                inPtr++;

            }
        }
        std::string B_BOX_STRING = std::to_string(bounding_box*2 + 1);
        imOut.writeFits(imageDestination + B_BOX_STRING + "x" + B_BOX_STRING + ".fits");
    #endif
}


//Kernels are not required to be normalized, but are not individually normalized.
//The total linear combination of all basis kernels is normalized.
//Individual kernel normalization can be controlled using the polynomial coefficients.
void generalKernel::createLinearCombinationProgram(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels){

    if(weights.size() != kernels.size()){
        cout << "ERROR, must supply equal number of weights and kernels" << endl;
        return;
    }

    //Condensed version generalized for any number of kernels/weights
    Expr blur_image_help = 0.0f;
    Expr blur_variance_help = 0.0f;
    Expr norm = 0.0f;
    Expr cur_kernel_location;
    total_mask_output = cast<uint16_t>(0);  //make sure blur_mask_help has type uint16

    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            cur_kernel_location = 0.0f;
            for(int h = 0; h < kernels.size(); h++){
                cur_kernel_location += weights[h]()*(*kernels[h])(i, j);
                norm += weights[h]()*(*kernels[h])(i, j);
            }
            blur_image_help += image_bounded(x + i, y + j) * cur_kernel_location;
            blur_variance_help += variance_bounded(x + i, y + j) 
                                    * cur_kernel_location * cur_kernel_location;

            total_mask_output = select(cur_kernel_location == 0.0f, total_mask_output,
                                total_mask_output | mask_bounded(x + i, y + j));
//            total_mask_output = total_mask_output | mask_bounded(x + i, y + j);  
        }
    }
    total_image_output = blur_image_help/norm;
    total_variance_output = blur_variance_help/(norm*norm);



    //Explicit version for 5 kernels/weights
    //Compute output image plane
/*    Expr blur_image_help = 0.0f;
    Expr norm = 0.0f;
    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            blur_image_help += image_bounded(x + i, y + j) * (weights[0]()*kernels[0](i, j) +
                weights[1]()*kernels[1](i, j) + weights[2]()*kernels[2](i, j) + 
                weights[3]()*kernels[3](i, j) + weights[4]()*kernels[4](i, j)); 
            norm += (weights[0]()*kernels[0](i, j) +
                weights[1]()*kernels[1](i, j) + weights[2]()*kernels[2](i, j) + 
                weights[3]()*kernels[3](i, j) + weights[4]()*kernels[4](i, j)); 
        }
    }
    total_image_output = blur_image_help/norm;

    //Compute output variance plane
    Expr blur_variance_help = 0.0f;
    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            blur_variance_help += variance_bounded(x + i, y + j) * (weights[0]()*kernels[0](i, j) +
                weights[1]()*kernels[1](i, j) + weights[2]()*kernels[2](i, j) + 
                weights[3]()*kernels[3](i, j) + weights[4]()*kernels[4](i, j))
                *(weights[0]()*kernels[0](i, j) +
                weights[1]()*kernels[1](i, j) + weights[2]()*kernels[2](i, j) + 
                weights[3]()*kernels[3](i, j) + weights[4]()*kernels[4](i, j)); 
        }
    }
    total_variance_output = blur_variance_help/(norm*norm);

    //Compute output mask plane
    total_mask_output = cast<uint16_t>(0);  //make sure blur_mask_help has type uint16
    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            total_mask_output = select((weights[0]()*kernels[0](i, j) + weights[1]()*kernels[1](i, j) +
            weights[2]()*kernels[2](i, j) + weights[3]()*kernels[3](i, j) +weights[4]()*kernels[4](i, j))
            == 0.0f, total_mask_output, total_mask_output | mask_bounded(x + i, y + j));
//            total_mask_output = total_mask_output | mask_bounded(x + i, y + j);    
        }
    }
*/
}

void generalKernelWithTuples::createLinearCombinationProgram(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels){
        generalKernel::createLinearCombinationProgram(weights, kernels);
        combined_output(x, y) = Tuple(total_image_output, total_variance_output, total_mask_output);
}



void generalKernelWithoutTuples::createLinearCombinationProgram(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels){
        generalKernel::createLinearCombinationProgram(weights, kernels);

        image_output(x, y) = total_image_output;
        variance_output(x, y) = total_variance_output;
        mask_output(x, y) = total_mask_output;
}


void generalKernelWithTuples::schedule_for_cpu() {
    // Split the y coordinate of the consumer into strips of 4 scanlines:
    combined_output.split(y, y_0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y_0);
    // Vectorize across x by a factor of four.
    combined_output.vectorize(x, 16);  
}

void generalKernelWithoutTuples::schedule_for_cpu() {
    image_output.split(y, y_0, yi, 32);
    image_output.parallel(y_0);
    image_output.vectorize(x, 16);

    variance_output.split(y, y_0, yi, 32);
    variance_output.parallel(y_0);
    variance_output.vectorize(x, 16);

    mask_output.split(y, y_0, yi, 32);
    mask_output.parallel(y_0);
    mask_output.vectorize(x, 16);  
}

void generalKernelWithTuples::schedule_for_gpu() {
    // Compute curved in 2D 8x8 tiles using the GPU.

    combined_output.gpu_tile(x, y, 16, 16);

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


void generalKernelWithoutTuples::schedule_for_gpu() {
    // Compute curved in 2D 8x8 tiles using the GPU.

    image_output.gpu_tile(x, y, 16, 16);
    variance_output.gpu_tile(x, y, 16, 16);
    mask_output.gpu_tile(x, y, 16, 16);

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

    image_output.compile_jit(target);
    variance_output.compile_jit(target);
    mask_output.compile_jit(target);
}

void generalKernelWithTuples::test_performance_cpu() {
    // Benchmark the pipeline.
    //Run once to initialize and save output
    Realization r = combined_output.realize(image.width(), image.height());
    image_output_cpu = r[0];
    variance_output_cpu = r[1];
    mask_output_cpu = r[2];

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
        t1 = current_time();
        r = combined_output.realize(image.width(), image.height());
        t2 = current_time();
        t3 = current_time();
        t4 = current_time();
        curTime = (t4-t1);
        
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


void generalKernelWithoutTuples::test_performance_cpu() {
    // Test the performance of the pipeline.


    //set the output dimensions equal to the input dimensions by copying input
    //planes into the output planes (image input is used for image and variance
    //output because the content doesn't matter, but the mask plane is a different type)
    image_output_cpu = image;
    variance_output_cpu = image;
    mask_output_cpu = mask;

    image_output.realize(image_output_cpu);
    variance_output.realize(variance_output_cpu);
    mask_output.realize(mask_output_cpu);


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

        t1 = current_time();
        image_output.realize(image_output_cpu);
        t2 = current_time();
        variance_output.realize(variance_output_cpu);
        t3 = current_time();
        mask_output.realize(mask_output_cpu);
        t4 = current_time();
        curTime = (t4-t1);
        
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

//test GPU performance and save 
void generalKernelWithTuples::test_performance_gpu() {
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

    // Run the filter once to initialize any GPU runtime state.
    // Also save the output image, variance, and mask planes
    Realization r = combined_output.realize(image.width(), image.height());
    image_output_gpu = r[0];
    variance_output_gpu = r[1];
    mask_output_gpu = r[2];

    // Now take the best of 3 runs for timing.
    double best_time;
    double t1;
    double t2;
    double elapsed;
    for (int i = 0; i < 3; i++) {

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

        elapsed = (t2 - t1)/100;
        if (i == 0 || elapsed < best_time) {
            best_time = elapsed;
        }
    }

    printf("%1.4f milliseconds\n", best_time);
}

void generalKernelWithoutTuples::test_performance_gpu() {
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

    image_output.realize(image_gpu_output);
    variance_output.realize(variance_gpu_output);
    mask_output.realize(mask_gpu_output);
    image_output_gpu = image_gpu_output;
    variance_output_gpu = variance_gpu_output;
    mask_output_gpu = mask_gpu_output;

    // Now take the best of 3 runs for timing.
    double best_time;
    double t1;
    double t2;
    double elapsed;
    for (int i = 0; i < 3; i++) {
        t1 = current_time();
        // Run the filter 100 times.
        for (int j = 0; j < 100; j++) {
            image_output.realize(image_gpu_output);
            variance_output.realize(variance_gpu_output);
            mask_output.realize(mask_gpu_output);
        }
        // Force any GPU code to finish by copying the buffer back to the CPU.
        image_gpu_output.copy_to_host();
        variance_gpu_output.copy_to_host();
        mask_gpu_output.copy_to_host();

        t2 = current_time();

        elapsed = (t2 - t1)/100;
        if (i == 0 || elapsed < best_time) {
            best_time = elapsed;
        }
    }

    printf("%1.4f milliseconds\n", best_time);
}



void generalKernel::debug(){
    //Check out what is happening
//    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
//    combined_output.compile_to_lowered_stmt("generalKernelBlurImage.html", {image}, HTML);
}

void generalKernel::test_correctness(string referenceLocation, int details) {
#ifndef STANDALONE //do nothing when STANDALONE is defined
    //load reference image, variance, and mask planes
    auto im = afwImage::MaskedImage<float>(referenceLocation);
    int width = im.getWidth(), height = im.getHeight();

    Image<float> reference_image(width, height);
    Image<float> reference_variance(width, height);
    Image<uint16_t> reference_mask(width, height);

    //Read image in
    for (int y = 0; y < height; y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = im.x_at(0, y);
        for (int x = 0; x < width; x++){
            reference_image(x, y) = (*inPtr).image();
            reference_variance(x, y) = (*inPtr).variance();
            reference_mask(x, y) = (*inPtr).mask();
            inPtr++;
        }
    }

    //compute the maximum differences of the input images' image, mask, and variance values
    //write an output image whose values are reference_image-image_output_cpu 
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
    int reference_maskVal = -1;
    int mask_output_cpuVal = -1;

    float maxImage1 = 0;
    float maxVariance1 = 0;
    int maxMask1 = 0;

    float maxImage2 = 0;
    float maxVariance2 = 0;
    int maxMask2 = 0;


    auto imDiff = afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>(width, height);
    for (int y = bounding_box; y < imDiff.getHeight()-bounding_box; y++) {
        afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>::x_iterator inPtr = imDiff.x_at(0, y);

        for (int x = bounding_box; x < imDiff.getWidth()-bounding_box; x++){
//            afwImage::pixel::SinglePixel<float, afwImage::MaskPixel, afwImage::VariancePixel> curPixel(reference_image(x, y) - image_output_cpu(x, y), reference_variance(x, y) - variance_output_cpu(x, y), reference_mask(x, y) - mask_output_cpu(x, y));
            afwImage::pixel::SinglePixel<float, afwImage::MaskPixel, afwImage::VariancePixel> curPixel((reference_image(x, y) - image_output_cpu(x, y))/reference_variance(x, y), reference_variance(x, y) - variance_output_cpu(x, y), reference_mask(x, y) - mask_output_cpu(x, y));
            (*inPtr) = curPixel;
            inPtr++;
//            if(abs(reference_image(x, y) - image_output_cpu(x, y))/min(abs(reference_image(x, y)), abs(image_output_cpu(x, y))) > maxImageDiffPercent){
            if((abs(reference_image(x, y) - image_output_cpu(x, y))/min(abs(reference_image(x, y)), abs(image_output_cpu(x, y)))) > maxImageDiff){
                maxImageDiffPercent = abs(reference_image(x, y) - image_output_cpu(x, y))/min(abs(reference_image(x, y)), abs(image_output_cpu(x, y)));
                maxImageDiff = abs(reference_image(x, y) - image_output_cpu(x, y));
                valueAtMaxImageDiff = min(abs(reference_image(x, y)), abs(image_output_cpu(x, y)));
                varianceAtMaxImgDiff = min(abs(reference_variance(x, y)), abs(variance_output_cpu(x, y)));
                imageX = x;
                imageY = y;
            }
            if(abs(reference_variance(x, y) - variance_output_cpu(x, y))/min(abs(reference_variance(x, y)), abs(variance_output_cpu(x, y))) > maxVarianceDiff){
                maxVarianceDiffPercent = abs(reference_variance(x, y) - variance_output_cpu(x, y))/min(abs(reference_variance(x, y)), abs(variance_output_cpu(x, y)));
                maxVarianceDiff = abs(reference_variance(x, y) - variance_output_cpu(x, y));
                valueAtMaxVarianceDiff1 = reference_variance(x, y);
                valueAtMaxVarianceDiff2 = variance_output_cpu(x, y);
                varianceX = x;
                varianceY = y;
            }
            if(abs(reference_mask(x, y) - mask_output_cpu(x, y)) > maxMaskDiff){
                maxMaskDiff = abs(reference_mask(x, y) - mask_output_cpu(x, y));
                maskX = x;
                maskY = y;
                reference_maskVal = reference_mask(x, y);
                mask_output_cpuVal = mask_output_cpu(x, y);
            }

            if(abs(reference_image(x, y)) > maxImage1)
                maxImage1 = abs(reference_image(x, y));
            if(abs(reference_variance(x, y)) > maxVariance1)
                maxVariance1 = abs(reference_variance(x, y));
            if(reference_mask(x, y) > maxMask1)
                maxMask1 = reference_mask(x, y);

            if(abs(image_output_cpu(x, y)) > maxImage2)
                maxImage2 = abs(image_output_cpu(x, y));
            if(abs(variance_output_cpu(x, y)) > maxVariance2)
                maxVariance2 = abs(variance_output_cpu(x, y));
            if(mask_output_cpu(x, y) > maxMask2)
                maxMask2 = mask_output_cpu(x, y);
        }
    }

    if(details == 0){
        cout << "Max (image difference)/(smaller of two image value) = " << maxImageDiffPercent << endl;
        cout << "Max (variance difference)/(smaller of two variance value) = " << maxVarianceDiffPercent << endl;
        cout << "Max mask difference = " << maxMaskDiff << endl;
    }
    else if(details ==1){
        cout << "Max (image difference)/(min img value) = " << maxImageDiffPercent << ",  Max image difference = " << maxImageDiff
        << ", value at max image difference = " << valueAtMaxImageDiff << " at position: (" << imageX << ", " << imageY
        << "), " << "variance at max difference = " << varianceAtMaxImgDiff << endl;

        cout << "Max (variance difference)/(min var value) = " << maxVarianceDiffPercent << ",  Max Variance difference = " << maxVarianceDiff
        << ", reference_variance at max Variance difference = " << valueAtMaxVarianceDiff1 << ", variance_output_cpu at max Variance difference = "
        << valueAtMaxVarianceDiff2 <<" at position: (" << varianceX << ", ";
        cout << varianceY << ")" << endl;

        cout << "Max mask difference = " << maxMaskDiff << " at position: (" << maskX << ", " << maskY << ")" << 
        ", img1 mask = " << reference_maskVal << ", img2 mask = " << mask_output_cpuVal << endl;

        cout << "Max reference_image = " << maxImage1 << ", max reference_variance = " << maxVariance1 << 
        ", max reference_mask = " << maxMask1 << endl;

        cout << "Max image_output_cpu = " << maxImage2 << ", max variance_output_cpu = " << maxVariance2 << 
        ", max mask_output_cpu = " << maxMask2 << endl; 
    }   
#endif
}

int main(int argc, char *argv[]) {


    //test a linear combination of 5 spatially invariant gaussians
    std::vector<polynomial> weights;
    for(int i = 1; i < 6; i++){
        polynomial curPol = polynomial(0.1f, 0.001f, 0.001f, 0.000001f, 0.000001f, 0.000001f,
                            0.000000001f, 0.000000001f, 0.000000001f, 0.000000001f);
//        polynomial curPol = polynomial(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
//                            0.0f, 0.0f, 0.0f, 0.0f);
        weights.push_back(curPol);
    }

    std::vector<kernel2D *> kernels;
    kernels.resize(5);
    gaussian2D k1 = gaussian2D(2.0f, 2.0f, 0.0f);
    gaussian2D k2 = gaussian2D(.5f, 4.0f, 0.0f);
    gaussian2D k3 = gaussian2D(.5f, 4.0f, M_PI/4.0f);
    gaussian2D k4 = gaussian2D(.5f, 4.0f, M_PI/2.0f);
    gaussian2D k5 = gaussian2D(4.0f, 4.0f, 0.0f);

    kernels[0] = &k1; 
    kernels[1] = &k2;
    kernels[2] = &k3;
    kernels[3] = &k4;
    kernels[4] = &k5;

    //create 5x5 kernel with 
    generalKernelWithTuples p1("./images/calexp-004207-g3-0123.fits", 5);
    p1.createLinearCombinationProgram(weights, kernels);
    p1.schedule_for_cpu();
    p1.test_performance_cpu();
    p1.test_correctness("./lsstLinearCombination5x5.fits", 0);

 //   cout << "Kernel size = " << (bounding_box*2 + 1) << " x " << (bounding_box*2 + 1) <<endl;

    //save fits image (does nothing if STANDALONE defined)
    //save_fits_image("folder1/folder2/imageName") will save image to:
    //folder1/folder2/imageNameKERNELSIZExKERNELSIZE.fits
    //e.g. for the case of a 5x5 kernel:
    //folder1/folder2/imageName5x5.fits
    p1.save_fits_image("./images/linearCombination/generalKernel");

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

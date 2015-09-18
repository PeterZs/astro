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
        printf("Kernel size %d x %d\n", kernelSize, kernelSize);
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


//Kernels are not required to be normalized, but are not individually normalized here.
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

            if(!OR_ALL_MASK_PIXELS){
                total_mask_output = select(cur_kernel_location == 0.0f, total_mask_output,
                                    total_mask_output | mask_bounded(x + i, y + j));
            }  
            else{ 
                total_mask_output = total_mask_output | mask_bounded(x + i, y + j);  
            }  
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


//Kernels are not required to be normalized, but are not individually normalized here.
//The total linear combination of all basis kernels is normalized at each point of the 
//original kernel evaluation, before interpolation.
//Individual kernel normalization can be controlled using the polynomial coefficients.
//2 dimensional linear interpolation is performed on the original kernel using a grid
//of size (interpDist x interpDist).
void generalKernel::createLinearCombinationProgramWithInterpolation(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels, int interpDist){

    if(weights.size() != kernels.size()){
        cout << "ERROR, must supply equal number of weights and kernels" << endl;
        return;
    }

    //Condensed version generalized for any number of kernels/weights
    Expr blur_image_help = 0.0f;
    Expr blur_variance_help = 0.0f;
    Expr norm = 0.0f;
    Expr cur_kernel_location = 0.0f;
    total_mask_output = cast<uint16_t>(0);  //make sure blur_mask_help has type uint16

    //The original linear combination kernel function, normalized after linear combination.
/*    Func normalizedTotalKernel;
    normalizedTotalKernel(x, y, i, j) = 0.0f;
    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            cur_kernel_location = 0.0f;
            for(int h = 0; h < kernels.size(); h++){
                cur_kernel_location += weights[h]()*(*kernels[h])(i, j);
                norm += weights[h]()*(*kernels[h])(i, j);
            }
            normalizedTotalKernel(x, y, i, j) = cur_kernel_location;
        }
    }
*/

    Func normalizedTotalKernel;

    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            for(int h = 0; h < kernels.size(); h++){
                norm += weights[h]()*(*kernels[h])(i, j);
            }
        }
    }

    Expr linear_combo_kernel = 0.0f;
    for(int h = 0; h < kernels.size(); h++){
        linear_combo_kernel += weights[h]()*(*kernels[h])();
    }
    normalizedTotalKernel(x, y, i, j) = linear_combo_kernel/norm;





    //This is a compressed version of the original kernel that can be evaluated
    //compute_root() without computing the kernel at points that will not
    //actually be used for interpolation.
    compressedKernel(x, y, i, j) = normalizedTotalKernel(x*interpDist, y*interpDist, i, j);


    //Do the actual interpolation now

//******************************Without the compressed kernel
/*    Func interpKernel;
    Expr xS = x%interpDist;
    Expr xG = interpDist - xS;
    Expr yS = y%interpDist;
    Expr yG = interpDist - yS;

    Expr xInterpS = ((xG * normalizedTotalKernel(x - xS, y - yS, i, j)) + (xS * normalizedTotalKernel(x + xG, y - yS, i, j)))/interpDist;
    Expr xInterpG = ((xG * normalizedTotalKernel(x - xS, y + yG, i, j)) + (xS * normalizedTotalKernel(x + xG, y + yG, i, j)))/interpDist;

    interpKernel(x, y, i, j) = ((yG * xInterpS) + (yS * xInterpG))/interpDist;
*/
//******************************DONE Without the compressed kernel


//******************************With the compressed kernel
    Func interpKernel;
    Expr xS = x%interpDist;
    Expr xG = interpDist - xS;
    Expr yS = y%interpDist;
    Expr yG = interpDist - yS;

    Expr xInterpS = ((xG * compressedKernel((x - xS)/interpDist, (y - yS)/interpDist, i, j)) +
                     (xS * compressedKernel((x + xG)/interpDist, (y - yS)/interpDist, i, j))) / interpDist;

    Expr xInterpG = ((xG * compressedKernel((x - xS)/interpDist, (y + yG)/interpDist, i, j)) +
                     (xS * compressedKernel((x + xG)/interpDist, (y + yG)/interpDist, i, j))) / interpDist;

    interpKernel(x, y, i, j) = ((yG * xInterpS) + (yS * xInterpG))/interpDist;
//******************************DONE With the compressed kernel

    //perform convolution with the interpolated kernel
    total_image_output = 0.0f;
    total_variance_output = 0.0f;
    total_mask_output = cast<uint16_t>(0);
    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            total_image_output += image_bounded(x + i, y + j)*interpKernel(x, y, i, j);
            total_variance_output += variance_bounded(x + i, y + j)*interpKernel(x, y, i, j)
                                    *interpKernel(x, y, i, j);

            if(!OR_ALL_MASK_PIXELS){
                total_mask_output = select(interpKernel(x, y, i, j) == 0.0f, total_mask_output,
                                    total_mask_output | mask_bounded(x + i, y + j));
            }  
            else{ 
                total_mask_output = total_mask_output | mask_bounded(x + i, y + j);  
            } 
        }
    }

}

//Kernels are not required to be normalized, but are not individually normalized here.
//The total linear combination of all basis kernels is normalized at each point of the 
//original kernel evaluation, before interpolation.
//Individual kernel normalization can be controlled using the polynomial coefficients.
//2 dimensional linear interpolation is performed on the original kernel using a grid
//of size (interpDist x interpDist).
void generalKernel::createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels, int interpDist){
    if(weights.size() != kernels.size()){
        cout << "ERROR, must supply equal number of weights and kernels" << endl;
        return;
    }

    Expr linear_combo_kernel = 0.0f;
    for(int h = 0; h < kernels.size(); h++){
        linear_combo_kernel += weights[h]()*(*kernels[h])();
    }

    Func kernel;
    kernel(x, y, i, j) = linear_combo_kernel;

    //This is a compressed version of the original kernel that can be evaluated
    //compute_root() without computing the kernel at points that will not
    //actually be used for interpolation.
    compressedKernel(x, y, i, j) = kernel(x*interpDist, y*interpDist, i, j);


    //Do the actual interpolation now

//******************************Without the compressed kernel
    Func interpKernel;
    Expr xS = x%interpDist;
    Expr xG = interpDist - xS;
    Expr yS = y%interpDist;
    Expr yG = interpDist - yS;

    Expr xInterpS = ((xG * kernel(x - xS, y - yS, i, j)) + (xS * kernel(x + xG, y - yS, i, j)))/interpDist;
    Expr xInterpG = ((xG * kernel(x - xS, y + yG, i, j)) + (xS * kernel(x + xG, y + yG, i, j)))/interpDist;

    interpKernel(x, y, i, j) = ((yG * xInterpS) + (yS * xInterpG))/interpDist;

//******************************DONE Without the compressed kernel


//******************************With the compressed kernel
/*    Func interpKernel;
    Expr xS = x%interpDist;
    Expr xG = interpDist - xS;
    Expr yS = y%interpDist;
    Expr yG = interpDist - yS;

    Expr xInterpS = ((xG * compressedKernel((x - xS)/interpDist, (y - yS)/interpDist, i, j)) +
                     (xS * compressedKernel((x + xG)/interpDist, (y - yS)/interpDist, i, j))) / interpDist;

    Expr xInterpG = ((xG * compressedKernel((x - xS)/interpDist, (y + yG)/interpDist, i, j)) +
                     (xS * compressedKernel((x + xG)/interpDist, (y + yG)/interpDist, i, j))) / interpDist;

    interpKernel(x, y, i, j) = ((yG * xInterpS) + (yS * xInterpG))/interpDist;
*/
//******************************DONE With the compressed kernel

    //perform convolution with the interpolated kernel
    total_image_output = 0.0f;
    total_variance_output = 0.0f;
    total_mask_output = cast<uint16_t>(0);
    Expr norm = 0.0f;
    for(int i = -bounding_box; i <= bounding_box; i++){
        for(int j = -bounding_box; j <= bounding_box; j++){
            total_image_output += image_bounded(x + i, y + j)*interpKernel(x, y, i, j);
            total_variance_output += variance_bounded(x + i, y + j)*interpKernel(x, y, i, j)
                                    *interpKernel(x, y, i, j);
            if(!OR_ALL_MASK_PIXELS){
                total_mask_output = select(interpKernel(x, y, i, j) == 0.0f, total_mask_output,
                                    total_mask_output | mask_bounded(x + i, y + j));
            }  
            else{ 
                total_mask_output = total_mask_output | mask_bounded(x + i, y + j);  
            } 
            norm += interpKernel(x, y, i, j);
        }
    }
    total_image_output = total_image_output/norm;
    total_variance_output = total_variance_output/(norm*norm);

}




void generalKernelWithTuples::createLinearCombinationProgramWithInterpolation(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels, int interpDist){
        generalKernel::createLinearCombinationProgramWithInterpolation(weights, kernels, interpDist);
        combined_output(x, y) = Tuple(total_image_output, total_variance_output, total_mask_output);
}

void generalKernelWithTuples::createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels, int interpDist){
        generalKernel::createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(weights, kernels, interpDist);
        combined_output(x, y) = Tuple(total_image_output, total_variance_output, total_mask_output);
}

void generalKernelWithoutTuples::createLinearCombinationProgramWithInterpolation(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels, int interpDist){
        generalKernel::createLinearCombinationProgramWithInterpolation(weights, kernels, interpDist);

        image_output(x, y) = total_image_output;
        variance_output(x, y) = total_variance_output;
        mask_output(x, y) = total_mask_output;
}

void generalKernelWithoutTuples::createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(std::vector<polynomial> weights,
    std::vector<kernel2D *> kernels, int interpDist){
        generalKernel::createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(weights, kernels, interpDist);

        image_output(x, y) = total_image_output;
        variance_output(x, y) = total_variance_output;
        mask_output(x, y) = total_mask_output;
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

void generalKernelWithTuples::schedule_interpolation_for_cpu() {
    // Split the y coordinate of the consumer into strips of 4 scanlines:
    combined_output.split(y, y_0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y_0);
    // Vectorize across x by a factor of four.
    combined_output.vectorize(x, 16);  

    compressedKernel.compute_root();
    // Split the y coordinate of the consumer into strips of 4 scanlines:
    compressedKernel.split(y, y_0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    compressedKernel.parallel(y_0);
    // Vectorize across x by a factor of four.
    compressedKernel.vectorize(x, 16);  
}

void generalKernelWithoutTuples::schedule_for_cpu() {
    image_output.split(y, y_0, yi, 10);
    image_output.parallel(y_0);
    image_output.vectorize(x, 8);

    variance_output.split(y, y_0, yi, 10);
    variance_output.parallel(y_0);
    variance_output.vectorize(x, 8);

    mask_output.split(y, y_0, yi, 10);
    mask_output.parallel(y_0);
    mask_output.vectorize(x, 8);  
}

void generalKernelWithoutTuples::schedule_interpolation_for_cpu() {
    image_output.split(y, y_0, yi, 10);
    image_output.parallel(y_0);
    image_output.vectorize(x, 8);

    variance_output.split(y, y_0, yi, 10);
    variance_output.parallel(y_0);
    variance_output.vectorize(x, 8);

    mask_output.split(y, y_0, yi, 10);
    mask_output.parallel(y_0);
    mask_output.vectorize(x, 8);  

    compressedKernel.compute_root();

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
    std::cout << "Using tuples, average Time: " << average << ", Min = " <<
    min << ", Max = " << max << ", with " << NUMBER_OF_RUNS <<
    " runs" << '\n';
}


void generalKernelWithoutTuples::test_performance_cpu() {
    // Test the performance of the pipeline.

//Maybe this is a problem?
/*    //set the output dimensions equal to the input dimensions by copying input
    //planes into the output planes (image input is used for image and variance
    //output because the content doesn't matter, but the mask plane is a different type)
    image_output_cpu = image;
    variance_output_cpu = image;
    mask_output_cpu = mask;
*/
    Image<float> allocateImageOut(image.width(), image.height());
    Image<float> allocateVarianceOut(image.width(), image.height());
    Image<uint16_t> allocateMaskOut(image.width(), image.height());
    image_output_cpu = allocateImageOut;
    variance_output_cpu = allocateVarianceOut;
    mask_output_cpu = allocateMaskOut;


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
    compressedKernel.trace_realizations();
//    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
//    combined_output.compile_to_lowered_stmt("generalKernelBlurImage.html", {image}, HTML);
}

void generalKernelWithTuples::debug(){
    //Check out what is happening
    compressedKernel.trace_realizations();
    combined_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
    combined_output.compile_to_lowered_stmt("generalKernelWithTuplesCombined_Output.html", {image}, HTML);
}

void generalKernelWithoutTuples::debug(){
    //Check out what is happening
    compressedKernel.trace_realizations();
    image_output.print_loop_nest();
    // Print out pseudocode for the pipeline.
    image_output.compile_to_lowered_stmt("generalKernelBlurImageWithoutTuplesImage_Output.html", {image}, HTML);
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


    int edgeDistance = bounding_box + 1;

    auto imDiff = afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>(width, height);
    for (int y = edgeDistance; y < imDiff.getHeight()-edgeDistance; y++) {
        afwImage::MaskedImage<float, afwImage::MaskPixel, afwImage::VariancePixel>::x_iterator inPtr = imDiff.x_at(0, y);

        for (int x = edgeDistance; x < imDiff.getWidth()-edgeDistance; x++){
//            afwImage::pixel::SinglePixel<float, afwImage::MaskPixel, afwImage::VariancePixel> curPixel(reference_image(x, y) - image_output_cpu(x, y), reference_variance(x, y) - variance_output_cpu(x, y), reference_mask(x, y) - mask_output_cpu(x, y));
            afwImage::pixel::SinglePixel<float, afwImage::MaskPixel, afwImage::VariancePixel> curPixel((reference_image(x, y) - image_output_cpu(x, y))/reference_variance(x, y), reference_variance(x, y) - variance_output_cpu(x, y), reference_mask(x, y) - mask_output_cpu(x, y));
            (*inPtr) = curPixel;
            inPtr++;
//            if(abs(reference_image(x, y) - image_output_cpu(x, y))/min(abs(reference_image(x, y)), abs(image_output_cpu(x, y))) > maxImageDiffPercent){
            if((abs(reference_image(x, y) - image_output_cpu(x, y))) > maxImageDiff){
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

            //look into mask in more detail
            if(details == 2){
                if(reference_mask(x, y) != mask_output_cpu(x, y)){
                    cout << "Halide and referene LSST disagree at (" << x << ", " << y << 
                    ")  Halide mask = " << mask_output_cpu(x, y) << ", LSST mask = " << reference_mask(x, y)
                    << endl;
                    cout <<"Original mask plane around problem pixel:" << endl;
                    for(int j = -bounding_box; j <=bounding_box; j++){
                        for(int i = -bounding_box; i <=bounding_box; i++){
                            cout << mask(x + i, y + j) << "\t";
                        }
                        cout << endl;
                    }
/*                    cout << endl << "Kernel around problem pixel: " << endl;
                    for(int j = -bounding_box; j <=bounding_box; j++){
                        for(int i = -bounding_box; i <=bounding_box; i++){
                            cout << kernelOut(x, y, i, j) << "\t";
                        }
                        cout << endl;
                    }
                    cout << endl << "cKernel around problem pixel: " << endl;
                    for(int j = -bounding_box; j <=bounding_box; j++){
                        for(int i = -bounding_box; i <=bounding_box; i++){
                            cout << cKernel(x, y, i, j) << "\t";
                        }
                        cout << endl;
                    }
*/
                }
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
        cout << "Max (image difference)/(smaller of two image values) = " << maxImageDiffPercent << endl;
        cout << "Max (variance difference)/(smaller of two variance values) = " << maxVarianceDiffPercent << endl;
        cout << "Max mask difference = " << maxMaskDiff << endl;
    }
    else if(details > 0){
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



void checkLinearComboAndAnalytic(){
    vector<int> kernelSizes;
    kernelSizes.resize(3);
    kernelSizes[0] = 5;
    kernelSizes[1] = 11;
    kernelSizes[2] = 19;

    for(int ii = 0; ii < kernelSizes.size(); ii++){
        //test a linear combination of 5 spatially invariant gaussians
        std::vector<polynomial> weights;
        for(int i = 1; i < 6; i++){
            polynomial curPol = polynomial(0.1f, 0.001f, 0.001f, 0.000001f, 0.000001f, 0.000001f,
                                0.000000001f, 0.000000001f, 0.000000001f, 0.000000001f);
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

        generalKernelWithTuples p1("./images/calexp-004207-g3-0123.fits", kernelSizes[ii]);
        p1.createLinearCombinationProgram(weights, kernels);
        p1.schedule_for_cpu();
        p1.test_performance_cpu();

        std::string curKernelSize = std::to_string(kernelSizes[ii]);
        p1.test_correctness("./lsstOutput/lsstLinearCombination" + curKernelSize + "x" + curKernelSize +
            ".fits", 0);

        //test a single analytic kernel
        cout << "testing analytic kernel" << endl;

        polynomial poly(0.1f, 0.0f, 0.0019476158495634653f, 0.000001f, 0.000001f, 0.000001f,
                        0.000000001f, 0.000000001f, 0.000000001f, 0.000000001f);
        spatiallyVaryingGaussian2D analyticKernel(poly, poly, poly);
        vector<polynomial> weights1;
        weights1.push_back(polynomial(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
        vector<kernel2D *> singleKernel;
        singleKernel.push_back(&analyticKernel);
       
        generalKernelWithTuples p2("./images/calexp-004207-g3-0123.fits", kernelSizes[ii]);
        p2.createLinearCombinationProgram(weights1, singleKernel);
        p2.schedule_for_cpu();
        p2.test_performance_cpu();
        p2.test_correctness("./lsstOutput/lsstAnalyticKernel" + curKernelSize + "x" + curKernelSize +
            ".fits", 0);
    }

}


void checkInterpolation(){
    //test a single analytic kernel

    polynomial poly(0.1f, 0.0f, 0.0019476158495634653f, 0.000001f, 0.000001f, 0.000001f,
                    0.000000001f, 0.000000001f, 0.000000001f, 0.000000001f);
    spatiallyVaryingGaussian2D analyticKernel(poly, poly, poly);
    vector<polynomial> weights1;
    weights1.push_back(polynomial(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
    vector<kernel2D *> singleKernel;
    singleKernel.push_back(&analyticKernel);
  
    cout << "testing linear interpolation of analytic kernel, normalization after interpolation, with tuples" << endl;

    generalKernelWithTuples p0("./images/calexp-004207-g3-0123.fits", 5); //create 5x5 kernel
    p0.createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(weights1, singleKernel, 10); //interpolate 10x10 grid
    p0.schedule_interpolation_for_cpu();
    p0.debug();
    p0.test_performance_cpu();

    std::string curKernelSize = std::to_string(5);
    p0.test_correctness("./lsstOutput/lsstAnalyticKernel" + curKernelSize + "x" + curKernelSize +
        ".fits", 1);
    cout << endl << endl;

    cout << "testing linear interpolation of analytic kernel, normalization before interpolation, with tuples" << endl;

    generalKernelWithTuples p1("./images/calexp-004207-g3-0123.fits", 5); //create 5x5 kernel
    p1.createLinearCombinationProgramWithInterpolation(weights1, singleKernel, 10);
    p1.schedule_interpolation_for_cpu();
    p1.debug();
    p1.test_performance_cpu();

    p1.test_correctness("./lsstOutput/lsstAnalyticKernel" + curKernelSize + "x" + curKernelSize +
        ".fits", 1);
    cout << endl << endl;




/*    cout << "testing linear interpolation of analytic kernel, normalization after interpolation, no tuples" << endl;
    generalKernelWithoutTuples p1("./images/calexp-004207-g3-0123.fits", 5); //create 5x5 kernel
    p1.createLinearCombinationProgramWithInterpolationNormalizeAfterInterpolation(weights1, singleKernel, 10); //interpolate 10x10 grid
    p1.schedule_interpolation_for_cpu();
    p1.debug();
    p1.test_performance_cpu();

    p1.test_correctness("./lsstOutput/lsstAnalyticKernel" + curKernelSize + "x" + curKernelSize +
        ".fits", 1);
    p1.save_fits_image("./gkInterpWithoutTuples");
*/

/*    cout << "testing linear interpolation of analytic kernel, normalization before interpolation" << endl;

    generalKernelWithTuples p2("./images/calexp-004207-g3-0123.fits", 5); //create 5x5 kernel
    p2.createLinearCombinationProgramWithInterpolation(weights1, singleKernel, 10); //interpolate 10x10 grid
    p2.schedule_interpolation_for_cpu();
    p2.debug();
    p2.test_performance_cpu();

    p2.test_correctness("./lsstOutput/lsstAnalyticKernel" + curKernelSize + "x" + curKernelSize +
        ".fits", 1);
*/
}


int main(int argc, char *argv[]) {
//    checkLinearComboAndAnalytic();
    checkInterpolation();

}


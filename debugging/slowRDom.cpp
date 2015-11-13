//on os x:
// g++ slowRDom.cpp -g -I ../include -L ../bin -lHalide -o slowRDom -std=c++11
// DYLD_LIBRARY_PATH=../bin ./slowRDom
//on linux:
// g++ slowRDom.cpp -g -I ../include -L ../bin -lHalide -lpthread -ldl -o slowRDom -std=c++11
// LD_LIBRARY_PATH=../bin ./slowRDom

#include <stdio.h>
#include "Halide.h"
#include <bitset>
#include "clock.h"
using namespace std;
using namespace Halide;

using Halide::Image;


int main(int argc, char *argv[]) {

    int width = 2048, height = 1489;
    printf("image is calculated without using an RDom, image2 is calculated usng an RDom\n");
    printf("Testing with dimensions: %d x %d\n", width, height);

    Image<float> image(width, height);
    Image<float> image2(width, height);

    int boundingBox = 2; 
    Var x, y, i_v, y0, yi;

    //compute output image and image2
    //Polynomials that define weights of spatially variant linear combination of 5 kernels
    Func polynomial1, polynomial2, polynomial3, polynomial4, polynomial5;
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

    //duplicates for RDom test
    Func polynomial1V, polynomial2V, polynomial3V, polynomial4V, polynomial5V;
    polynomial1V(x, y) = 0.1f + 0.002f*x + 0.003f*y + 0.4f*x*x + 0.5f*x*y
                     + 0.6f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
                     + 0.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial2V(x, y) = 1.1f + 1.002f*x + 1.003f*y + 1.4f*x*x + 1.5f*x*y
                     + 1.6f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
                     + 1.00011f*y*y*y;

    //for experimenting with optimizations

    polynomial3V(x, y) = 2.1f + 2.002f*x + 2.003f*y + 2.4f*x*x + 2.5f*x*y
                     + 2.6f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
                     + 2.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial4V(x, y) = 3.1f + 3.002f*x + 3.003f*y + 3.4f*x*x + 3.5f*x*y
                     + 3.6f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
                     + 3.00011f*y*y*y;

    //for experimenting with optimizations
    polynomial5V(x, y) = 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
                     + 4.00011f*y*y*y;

    //Kernel #1
    Func kernel1;
    float sigmaX1 = 2.0f;
    float sigmaY1 = 2.0f;
    float theta1 = 0.0f; //rotation of sigmaX axis
    kernel1(x, y) = (exp(-((x*cos(theta1) +y*sin(theta1))*(x*cos(theta1) +y*sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (sqrtf(2*M_PI)*sigmaX1))
                    *(exp(-((y*cos(theta1) - x*sin(theta1))*(y*cos(theta1) - x*sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (sqrtf(2*M_PI)*sigmaY1));


    //Kernel #2
    Func kernel2;
    float sigmaX2 = 0.5f;
    float sigmaY2 = 4.0f;
    float theta2 = 0.0f; //rotation of sigmaX axis
    kernel2(x, y) = (exp(-((x*cos(theta2) +y*sin(theta2))*(x*cos(theta2) +y*sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (sqrtf(2*M_PI)*sigmaX2))
                    *(exp(-((y*cos(theta2) - x*sin(theta2))*(y*cos(theta2) - x*sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (sqrtf(2*M_PI)*sigmaY2));

    //Kernel #3
    Func kernel3;
    float sigmaX3 = 0.5f;
    float sigmaY3 = 4.0f;
    float theta3 = 3.14159f/4; //rotation of sigmaX axis
    kernel3(x, y) = (exp(-((x*cos(theta3) +y*sin(theta3))*(x*cos(theta3) +y*sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (sqrtf(2*M_PI)*sigmaX3))
                    *(exp(-((y*cos(theta3) - x*sin(theta3))*(y*cos(theta3) - x*sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (sqrtf(2*M_PI)*sigmaY3));
    //Kernel #4
    Func kernel4;
    float sigmaX4 = 0.5f;
    float sigmaY4 = 4.0f;
    float theta4 = 3.14159f/2; //rotation of sigmaX axis
    kernel4(x, y) = (exp(-((x*cos(theta4) +y*sin(theta4))*(x*cos(theta4) +y*sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (sqrtf(2*M_PI)*sigmaX4))
                    *(exp(-((y*cos(theta4) - x*sin(theta4))*(y*cos(theta4) - x*sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (sqrtf(2*M_PI)*sigmaY4));


    //Kernel #5
    Func kernel5;
    float sigmaX5 = 4.0f;
    float sigmaY5 = 4.0f;
    float theta5 = 0.0; //rotation of sigmaX axis
    kernel5(x, y) = (exp(-((x*cos(theta5) +y*sin(theta5))*(x*cos(theta5) +y*sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (sqrtf(2*M_PI)*sigmaX5))
                    *(exp(-((y*cos(theta5) - x*sin(theta5))*(y*cos(theta5) - x*sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (sqrtf(2*M_PI)*sigmaY5));


    //Compute output image plane
    Func image_bounded ("image_bounded");
    image_bounded = BoundaryConditions::repeat_edge(image);
    Expr blur_image_help = 0.0f;
    Expr norm = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            blur_image_help += image_bounded(x + i, y + j) * (polynomial1(x, y)*kernel1(i, j) +
                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j)); 
            norm += (polynomial1(x, y)*kernel1(i, j) + polynomial2(x, y)*kernel2(i, j) + 
                polynomial3(x, y)*kernel3(i, j) + polynomial4(x, y)*kernel4(i, j) + 
                polynomial5(x, y)*kernel5(i, j));
        }
    }
    blur_image_help = blur_image_help/norm;
    Func blurImage ("blurImage");
    blurImage(x, y) = blur_image_help;





    //Compute output image2 plane
    Func image2_bounded ("image2_bounded");
    image2_bounded = BoundaryConditions::repeat_edge(image2);
    Func blurimage2 ("blurimage2");
    Func norm2 ("norm2");
    Expr blur_image2_help = 0.0f;

    RDom r(-boundingBox, 2*boundingBox+1, -boundingBox, 2*boundingBox+1);

    norm2(x, y) = sum(polynomial1V(x, y)*kernel1(r.x, r.y) +
                polynomial2V(x, y)*kernel2(r.x, r.y) + polynomial3V(x, y)*kernel3(r.x, r.y) + 
                polynomial4V(x, y)*kernel4(r.x, r.y) + polynomial5V(x, y)*kernel5(r.x, r.y));

    blurimage2(x, y) = sum(image2_bounded(x + r.x, y + r.y) * (polynomial1V(x, y)*kernel1(r.x, r.y) +
                polynomial2V(x, y)*kernel2(r.x, r.y) + polynomial3V(x, y)*kernel3(r.x, r.y) + 
                polynomial4V(x, y)*kernel4(r.x, r.y) + polynomial5V(x, y)*kernel5(r.x, r.y)));

    blurimage2(x, y) = blurimage2(x, y)/norm2(x, y);
//    polynomial1V.compute_at(blurimage2, x);
//    polynomial2V.compute_at(blurimage2, x);
//    polynomial3V.compute_at(blurimage2, x);
//    polynomial4V.compute_at(blurimage2, x);
//    polynomial5V.compute_at(blurimage2, x);

    polynomial1V.compute_root();
    polynomial2V.compute_root();
    polynomial3V.compute_root();
    polynomial4V.compute_root();
    polynomial5V.compute_root();

    kernel1.compute_root();
    kernel2.compute_root();
    kernel3.compute_root();
    kernel4.compute_root();
    kernel5.compute_root();


    //Check out what is happening
//    blurImage.print_loop_nest();
    // Print out pseudocode for the pipeline.
    blurImage.compile_to_lowered_stmt("slowRDomImage.html", {image}, HTML);
    blurImage.compile_to_c("slowRDomImage_C_Code.cpp", std::vector<Argument>(), "slowRDomImage_C_Code");

    blurimage2.compile_to_lowered_stmt("slowRDomImage2.html", {image}, HTML);
    blurimage2.compile_to_c("slowRDomImage2_C_Code.cpp", std::vector<Argument>(), "slowRDomImage2_C_Code");
//    blurimage2.compile_to_lowered_stmt("blur.html", {image2}, HTML);




    Image<float> image_output(image.width(), image.height());

    blurImage.realize(image_output);

    Image<float> image2_output(image2.width(), image2.height());
    blurimage2.realize(image2_output);


	double average = 0;
    double min;
    double max;
    double imgTime;
    double varTime;
    double maskTime;
    int numberOfRuns = 1;

    double minImg;
    double minVar;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        blurImage.realize(image_output);
        double t2 = current_time();
        blurimage2.realize(image2_output);
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
            minImg = imgTime;
            minVar = varTime;
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

            if(t2-t1 < minImg)
                minImg = t2-t1;
            if(t3-t2 < minVar)
                minVar = t3-t2;
        }
    }
    average = average/numberOfRuns;
//    std::cout << "Average Time: " << average << ", Min = " <<
//    min << ", Max = " << max << ", with " << numberOfRuns <<
//    " runs" << '\n';
//    cout << "For fastest run total time = " << min << ", imgTime = " << imgTime << ", varTime = " << varTime << 
//    "maskTime = " << maskTime << endl;
    printf("image time = %f, image2 time = %f\n", minImg, minVar);


}



#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>

#include <cstdio>
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>

#include <png.h>

class RandomVectorGenerator
{
private:
    std::mt19937_64 generator;
    std::normal_distribution<> dist;

public:
    RandomVectorGenerator() {}
    RandomVectorGenerator(std::mt19937_64::result_type seed): generator(seed) {}

    glm::dvec3 operator()() noexcept
    {
        glm::dvec3 result;
        result.x = dist(generator);
        result.y = dist(generator);
        result.z = dist(generator);
        return glm::normalize(result);
    }
};

template <typename T>
class Grid3D
{
private:
    std::size_t dimensions[3];
    T *data;

public:
    Grid3D(std::size_t x, std::size_t y, std::size_t z)
    {
        // Let's hope there's no overflows...
        data = new T[z * y * x + y * x + x];
        dimensions[0] = x;
        dimensions[1] = y;
        dimensions[2] = z;
    }
    Grid3D(const Grid3D& other)
    {
        std::size_t s = 0;
        std::size_t sAccum = 1;
        for(unsigned int i = 0; i < 3; ++i)
        {
            dimensions[i] = other.dimensions[i];
            sAccum *= dimensions[i];
            s += sAccum;
        }
        data = new T[s];
    }
    Grid3D(Grid3D&& other)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            dimensions[i] = other.dimensions[i];
            other.dimensions[i] = 0;
        }
        data = other.data;
        other.data = nullptr;
    }
    ~Grid3D()
    {
        delete[] data;
    }

    std::size_t getX() const
    {
        return dimensions[0];
    }
    std::size_t getY() const
    {
        return dimensions[1];
    }
    std::size_t getZ() const
    {
        return dimensions[2];
    }

    const T& index(std::size_t x, std::size_t y, std::size_t z) const
    {
        if(x >= dimensions[0])
            x = 0;
        if(y >= dimensions[1])
            y = 0;
        if(z >= dimensions[2])
            z = 0;
        return data[z * dimensions[1] * dimensions[0] + y * dimensions[0] + x];
    }
    T& index(std::size_t x, std::size_t y, std::size_t z)
    {
        if(x >= dimensions[0])
            x = 0;
        if(y >= dimensions[1])
            y = 0;
        if(z >= dimensions[2])
            z = 0;
        return data[z * dimensions[1] * dimensions[0] + y * dimensions[0] + x];
    }

    const T& indexRelative(std::size_t x, std::size_t y, std::size_t z,
                           std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) const
    {
        if(i >= dimensions[0] - x)
            i = dimensions[0] - x - 1;
        else if(-i >= x)
            i = x = 0;

        if(j >= dimensions[1] - y)
            j = dimensions[1] - y - 1;
        else if(-j >= y)
            j = y = 0;

        if(k >= dimensions[2] - z)
            k = dimensions[2] - z - 1;
        else if(-k >= z)
            k = z = 0;

        return data[k * y * x + j * x + i];
    }
    T& indexRelative(std::size_t x, std::size_t y, std::size_t z,
                     std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k)
    {
        return const_cast<T&>(indexRelative(x, y, z, i, j, k));
    }
};

const double pi = 3.14159265d;

template <class Grid>
class Perlin3D
{
private:
    const Grid& grid;

public:
    Perlin3D(const Grid& grid_in): grid(grid_in) {}

    double sample(const glm::dvec3& point) const
    {
        // x, y, z of index into grid for the cube that contains the point
        // Found by simply truncating (flooring since this is always positive)
        // the point coordinates
        std::size_t cubeOrigin[3] = {point.x, point.y, point.z};

        auto sample =
        [&cubeOrigin, &point, this](std::size_t x, std::size_t y, std::size_t z) -> double
        {
            glm::dvec3 grad = grid.index(cubeOrigin[0]+x, cubeOrigin[1]+y, cubeOrigin[2]+z);
            glm::dvec3 dist = point - glm::dvec3(cubeOrigin[0]+x, cubeOrigin[1]+y, cubeOrigin[2]+z);
            return glm::dot(dist, grad);
        };

        glm::dvec3 originDistance(point.x - cubeOrigin[0],
                                  point.y - cubeOrigin[1],
                                  point.z - cubeOrigin[2]);

        double xAverage1 = sample(0, 0, 0) * (1.0 - originDistance.x) +
                           sample(1, 0, 0) * originDistance.x;
        double xAverage2 = sample(0, 1, 0) * (1.0 - originDistance.x) +
                           sample(1, 1, 0) * originDistance.x;
        double xAverage3 = sample(0, 0, 1) * (1.0 - originDistance.x) +
                           sample(1, 0, 1) * originDistance.x;
        double xAverage4 = sample(0, 1, 1) * (1.0 - originDistance.x) +
                           sample(1, 1, 1) * originDistance.x;

        double yAverage1 = xAverage1 * (1.0 - originDistance.y) +
                           xAverage2 * originDistance.y;
        double yAverage2 = xAverage3 * (1.0 - originDistance.y) +
                           xAverage4 * originDistance.y;

        return yAverage1 * (1.0 - originDistance.z) + yAverage2 * originDistance.z;

//        double result = 0;
//        glm::dvec3 tempDistance;
//        for(std::size_t z = cubeOrigin[2]; z < cubeOrigin[2] + 2; ++z)
//        {
//            tempDistance.z = (double) z - point.z;
//            double yAverage = 0;
//            for(std::size_t y = cubeOrigin[1]; y < cubeOrigin[1] + 2; ++y)
//            {
//                tempDistance.y = (double) y - point.y;
//                double xAverage = 0;
//                for(std::size_t x = cubeOrigin[0]; x < cubeOrigin[0] + 2; ++x)
//                {
//                    tempDistance.x = (double) x - point.x;
//
//                    double dot = glm::dot(grid.index(x, y, z), tempDistance);
//                    // Take weighted average of dot products for x and x+1
//                    xAverage += dot * (1.0d - abs(tempDistance.x));
////                    xAverage += dot * pow(sin(abs(tempDistance.x)*pi/2.0d), 2);
//                }
//                // Take weighted averages of x averages for y and y+1
//                yAverage += xAverage * (1.0d - abs(tempDistance.y));
////                yAverage += xAverage * pow(sin(abs(tempDistance.y)*pi/2.0d), 2);
//            }
//            // Final result is weighted average of y averages for z and z+1
//            result += yAverage * (1.0d - abs(tempDistance.z));
////            result += yAverage * pow(sin(abs(tempDistance.z)*pi/2.0d), 2);
//        }
    }
};

int main()
{
    using namespace std;

    ofstream logger("log.txt");
    FILE *image = fopen("test.png", "wb");
    if(!image)
        return 1;

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                  nullptr, nullptr, nullptr);
    if(!png_ptr)
        return 1;

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr)
    {
        png_destroy_write_struct(&png_ptr, nullptr);
        return 1;
    }

    if(setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_write_struct(&png_ptr, nullptr);
        return 1;
    }

    png_init_io(png_ptr, image);
    png_set_IHDR(png_ptr, info_ptr, 4096, 4096, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    RandomVectorGenerator vecGen;

    Grid3D<glm::dvec3> test(64, 64, 64);
    for(unsigned int z = 0; z < 64; ++z)
    {
        for(unsigned int y = 0; y < 64; ++y)
        {
            for(unsigned int x = 0; x < 64; ++x)
            {
                test.index(x, y, z) = vecGen();
                logger << x << ' ' << y << ' ' << z << ": " << test.index(x, y, z) << endl;
            }
        }
    }

    Perlin3D<Grid3D<glm::dvec3>> perlin(test);
    vector<unsigned char> imageRows(4096 * 4096);
    vector<unsigned char*> imageRowPtrs(4096);
    for(unsigned int i = 0; i < 4096; ++i)
    {
        imageRowPtrs[i] = &(imageRows[i*4096]);
    }
    double largest = 0;
    double smallest = 1000;
    for(unsigned int y = 0; y < 4096; ++y)
    {
        for(unsigned int x = 0; x < 4096; ++x)
        {
            glm::dvec3 point(((double) x / 4096.d) * 64.0d, ((double) y / 4096.d) * 64.0d, 1.1d);
            double sample = 135.0d * abs(perlin.sample(point) + 0.5);
            imageRows[y * 4096 + x] = sample;
            largest = sample > largest ? sample : largest;
            smallest = sample < smallest ? sample : smallest;
        }
    }
    png_set_rows(png_ptr, info_ptr, &(imageRowPtrs[0]));
    png_write_png(png_ptr, info_ptr, 0, nullptr);
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    cout << largest << ' ' << smallest << endl;

//    FILE *dump = fopen("C:\\Development\\dump.txt", "w");
//    char *fileBuffer = new char[67108864];
//    setvbuf(dump, fileBuffer, _IOFBF, 67108864);
//
//    for(unsigned int z = 0; z < 512; ++z)
//    {
//        for(unsigned int y = 0; y < 512; ++y)
//        {
//            for(unsigned int x = 0; x < 512; ++x)
//            {
//                glm::dvec3& vec = test.index(x, y, z);
//                fprintf(dump, "(%f,%f,%f)\n", vec.x, vec.y, vec.z);
//            }
//        }
//    }
//
//    fclose(dump);
//    delete[] fileBuffer;
}

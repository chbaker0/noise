#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <limits>
#include <map>
#include <tuple>

#include <cstdio>
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtx/color_space.hpp>

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include <xxhash.h>

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

    void distanceTransform(glm::dvec3& originDistance) const
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            originDistance[i] = 3.0 * originDistance[i] * originDistance[i] -
                                2.0 * originDistance[i] * originDistance[i] * originDistance[i];
//            originDistance[i] = glm::cos(pi * originDistance[i]) / -2 + 0.5d;
        }
    }

public:
    Perlin3D(const Grid& grid_in): grid(grid_in) {}

    void getSamples(const glm::dvec3& start, const glm::dvec3& end, const glm::u64vec3& counts, double *buffer)
    {
        glm::u64vec3 currentSampleNum(0, 0, 0);
        glm::u64vec3 currentCube(0, 0, 0);

        glm::dvec3 sampleOffset;
        for(unsigned int i = 0; i < 2; ++i)
        {
            sampleOffset[i] = (start[i] - end[i]) / (double) counts[i];
        }

        glm::dvec3 grads[8];

        auto getGrads =
        [&grads, &currentCube, this]() -> void
        {
            for(unsigned int z = 0; z < 2; ++z)
                for(unsigned int y = 0; y < 2; ++y)
                    for(unsigned int x = 0; x < 2; ++x)
                        grads[4 * z + 2 * y + x] = grid.index(x, y, z);
        };

        getGrads();

        auto sample =
        [&grads, this](const glm::dvec3& point, std::size_t x, std::size_t y, std::size_t z) -> double
        {
            const glm::dvec3& grad = grads[4 * z + 2 * y + x];
            glm::dvec3 dist = point - glm::dvec3(x, y, z);
            return glm::dot(dist, grad);
        };

        do
        {
            glm::dvec3 samplePos;
            for(unsigned int i = 0; i < 2; ++i)
                samplePos[i] = start[i] + (double) currentSampleNum[i] * sampleOffset[i];

            glm::dvec3 originDistance(samplePos.x - currentCube[0],
                                      samplePos.y - currentCube[1],
                                      samplePos.z - currentCube[2]);
            distanceTransform(originDistance);

            glm::u64vec3 nextCube(samplePos.x, samplePos.y, samplePos.z);

            if(nextCube != currentCube)
            {
                currentCube = nextCube;
                getGrads();
            }

            double xAverage1 = sample(samplePos, 0, 0, 0) * (1.0 - originDistance.x) +
                               sample(samplePos, 1, 0, 0) * originDistance.x;
            double xAverage2 = sample(samplePos, 0, 1, 0) * (1.0 - originDistance.x) +
                               sample(samplePos, 1, 1, 0) * originDistance.x;
            double xAverage3 = sample(samplePos, 0, 0, 1) * (1.0 - originDistance.x) +
                               sample(samplePos, 1, 0, 1) * originDistance.x;
            double xAverage4 = sample(samplePos, 0, 1, 1) * (1.0 - originDistance.x) +
                               sample(samplePos, 1, 1, 1) * originDistance.x;

            double yAverage1 = xAverage1 * (1.0 - originDistance.y) +
                               xAverage2 * originDistance.y;
            double yAverage2 = xAverage3 * (1.0 - originDistance.y) +
                               xAverage4 * originDistance.y;

            *buffer = yAverage1 * (1.0 - originDistance.z) + yAverage2 * originDistance.z;

            currentSampleNum.x++;
            if(currentSampleNum.x >= counts.x)
                currentSampleNum.x = 0, currentSampleNum.y++;
            if(currentSampleNum.y >= counts.y)
                currentSampleNum.y = 0, currentSampleNum.z++;
        } while(currentSampleNum.z < counts.z);
    }

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
        distanceTransform(originDistance);

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
    }
};

class Grid3DComp
{
private:
    std::vector<std::mt19937_64::result_type> seedBuffer;

public:
    Grid3DComp(std::size_t s = 512)
    {
        seedBuffer.resize(s);

        std::mt19937_64 randGen;
        for(unsigned int i = 0; i < 512; ++i)
            randGen();
        for(std::size_t i = 0; i < seedBuffer.size(); ++i)
            seedBuffer[i] = randGen();
    }

    glm::dvec3 index(std::size_t x, std::size_t y, std::size_t z) const
    {
        std::size_t buffer[3] = {x, y, z};

        glm::dvec3 result;
        unsigned int i = 0;
        do
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                unsigned long long seed = seedBuffer[((x + 3) * (y + 7) + z + i*j) % seedBuffer.size()];
                unsigned long long hash = XXH64(buffer, sizeof(buffer), seed);

                using lim = std::numeric_limits<unsigned long long>;
                result[j] = (double) hash / (double) lim::max();
            }
        } while(i++ < 8 && glm::length(result) > 1.0d);

        return glm::normalize(result);
    }
};

int main()
{
    using namespace std;

    sf::RenderWindow win(sf::VideoMode(512, 512), "Perlin");
    if(!win.isOpen())
        return 1;

    RandomVectorGenerator vecGen;

    Grid3D<glm::dvec3> test(4, 4, 4);
    Grid3D<glm::dvec3> colors(4, 4, 4);
    for(unsigned int z = 0; z < 4; ++z)
    {
        for(unsigned int y = 0; y < 4; ++y)
        {
            for(unsigned int x = 0; x < 4; ++x)
            {
                test.index(x, y, z) = vecGen();
                colors.index(x, y, z) = vecGen();
            }
        }
    }

    Grid3DComp gridTest;

    Perlin3D<Grid3D<glm::dvec3>> perlin(test);
    Perlin3D<Grid3D<glm::dvec3>> perlinColor(colors);

    sf::Texture perlinTexture;

    if(!perlinTexture.create(512, 512))
        return 1;
    sf::Uint8 *pixels = new sf::Uint8[512 * 512 * 4];
    double *sampleBuffer = new double[512 * 512];

    sf::Sprite perlinSprite(perlinTexture);

    size_t frameCounter = 0;
    char filename[128];
    string filenameStr;
    while(win.isOpen())
    {
        sf::Event e;
        while(win.pollEvent(e))
        {
            if(e.type == sf::Event::Closed)
                win.close();
        }

        win.clear(sf::Color::Black);

//        perlin.getSamples(glm::dvec3(0.0, 0.0, frameCounter / 300.0d), glm::dvec3(4.0, 4.0, frameCounter / 300.0d), glm::u64vec3(512, 512, 1), sampleBuffer);

        #pragma omp parallel for
        for(unsigned int y = 0; y < 512; ++y)
        {
            for(unsigned int x = 0; x < 512; ++x)
            {
//                pixels[y * 512 * 4 + x * 4] = sampleBuffer[512 * y + x] * 127.0 + 127.0;
//                pixels[y * 512 * 4 + x * 4 + 1] = sampleBuffer[512 * y + x] * 127.0 + 127.0;
//                pixels[y * 512 * 4 + x * 4 + 2] = sampleBuffer[512 * y + x] * 127.0 + 127.0;
//                pixels[y * 512 * 4 + x * 4 + 3] = 255;
                glm::dvec3 point(((double) x / 512.d) * 4.0d, ((double) y / 512.d) * 4.0d, fmod(frameCounter / 300.0d, 4.0d));
                double sample = abs(perlin.sample(point) + 0.5);
                double colorSample = abs(perlinColor.sample(point) + 0.5);
                glm::dvec3 color = glm::rgbColor(glm::dvec3(colorSample * 360.0d, 1.0d, sample));
//                glm::dvec3 color(sample, sample, sample);
                pixels[y * 512 * 4 + x * 4] = glm::clamp(color.r * 255.0, 0.0, 254.9);
                pixels[y * 512 * 4 + x * 4 + 1] = glm::clamp(color.g * 255.0, 0.0, 254.9);
                pixels[y * 512 * 4 + x * 4 + 2] = glm::clamp(color.b * 255.0, 0.0, 254.9);
                pixels[y * 512 * 4 + x * 4 + 3] = 255;
            }
        }

        perlinTexture.update(pixels);

        perlinSprite.setPosition(sf::Vector2f(0, 0));
        win.draw(perlinSprite);

        win.display();

        ++frameCounter;
    }
}

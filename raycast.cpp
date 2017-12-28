#include <algorithm>
#include <cmath>

// Okay, it's ugly, I get it.
//
// I didn't use classes or a vector library.  I use loops to
// 3 everywhere for vector operations, and I even write out dot
// products.  My intent was to try to build a minimal triangle
// raycaster that I could convert simply to BASIC or FORTH or other
// simple languages without having to downconvert C++ concepts.
//
// I'm already thinking function calls are problematic, as well as
// the arrays with 3 dimensions.  I might have to inline the functions
// and use goto instead of early returns.
//
// So it's probably still not ugly enough.

const float F = 1000000;

const int W = 256;
const int H = 256;

const int M = 400000;

float tv[M][3][3]; // vertices
float tp[M][4]; // planes
float te[M][3][4]; // edge half-space planes
float tb[M][4]; // xmin, xmax, ymin, ymax

int tc = 0;

float L[3] = {.577, .577, .577};

void prepare_triangle(int k)
{
    // Get edge vectors
    float e[3][3];
    for(int i = 0; i < 3; i++)
        e[0][i] = tv[k][1][i] - tv[k][0][i];
    for(int i = 0; i < 3; i++)
        e[1][i] = tv[k][2][i] - tv[k][1][i];
    for(int i = 0; i < 3; i++)
        e[2][i] = tv[k][0][i] - tv[k][2][i];

    // Cross product e0 and e1 to get normal
    tp[k][0] = e[0][1] * e[1][2] - e[0][2] * e[1][1];
    tp[k][1] = e[0][2] * e[1][0] - e[0][0] * e[1][2];
    tp[k][2] = e[0][0] * e[1][1] - e[0][1] * e[1][0];

    // Normalize normal
    float d = sqrt(tp[k][0] * tp[k][0] + tp[k][1] * tp[k][1] + tp[k][2] * tp[k][2]);
    for(int i = 0; i < 3; i++)
        tp[k][i] /= d;

    // Use dot product of normal and v0 to get fourth plane equation element
    tp[k][3] = tp[k][0] * tv[k][0][0] + 
        tp[k][1] * tv[k][0][1] + 
        tp[k][2] * tv[k][0][2];

    // Get half-spaces representing each edge by making vector
    // perpendicular to edge and dot that with opposite vertex to get
    // plane equation
    for(int j = 0; j < 3; j++) {
        te[k][j][0] = tp[k][1] * e[j][2] - tp[k][2] * e[j][1];
        te[k][j][1] = tp[k][2] * e[j][0] - tp[k][0] * e[j][2];
        te[k][j][2] = tp[k][0] * e[j][1] - tp[k][1] * e[j][0];
        te[k][j][3] =
            te[k][j][0] * tv[k][j][0] +
            te[k][j][1] * tv[k][j][1] +
            te[k][j][2] * tv[k][j][2];
    }

    tb[k][0] = 2;
    tb[k][1] = -2;
    tb[k][2] = 2;
    tb[k][3] = -2;
    for(int j = 0; j < 3; j++) {
        tb[k][0] = std::min(tb[k][0], tv[k][j][0] / -tv[k][j][2]);
        tb[k][1] = std::max(tb[k][1], tv[k][j][0] / -tv[k][j][2]);
        tb[k][2] = std::min(tb[k][2], tv[k][j][1] / -tv[k][j][2]);
        tb[k][3] = std::max(tb[k][3], tv[k][j][1] / -tv[k][j][2]);
    }
}

float intersect_triangle(int k, float ray[3])
{
    int i;

    // find distance from ray origin (0,0,0) to plane in units of plane normal
    float factor = ray[0] * tp[k][0] + 
        ray[1] * tp[k][1] + 
        ray[2] * tp[k][2];

    // if coin, fail
    if(factor == 0.0)
       return F;

    // find distance to plane in world units
    float t = tp[k][3] / factor;

    // if intersection is behind the ray, fail
    if(t < 0)
        return F;

    // calculate the point in the plane
    float point[3];
    for(int i = 0; i < 3; i++)
        point[i] = ray[i] * t;

    for(i = 0; i < 3; i++) {
        // calculate the distance from each edge i to the point
        float edge_dot =
            te[k][i][0] * point[0] + 
            te[k][i][1] * point[1] + 
            te[k][i][2] * point[2];

        // if the distance is greater than the distance to the opposite vertex, fail
        if(edge_dot < te[k][i][3])
	   return F;
    }

    // succeed by returning the distance
    return t;
}

void reflect(float in[3], float n[3], float refl[3])
{
    float projection = in[0] * n[0] + in[1] * n[1] + in[2] * n[2];

    for(int i = 0; i < 3; i++)
        refl[i] = in[i] - 2.0f * projection * n[i];
}

int main()
{
    float v[3][3];

    // Read triangles.
    while(scanf("%f %f %f %f %f %f %f %f %f",
        &v[0][0], &v[0][1], &v[0][2],
        &v[1][0], &v[1][1], &v[1][2],
        &v[2][0], &v[2][1], &v[2][2]) == 9) {

        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++) {
                tv[tc][j][i] = v[j][i];
            }

        tc++;
    }

    // Get the bounding box of the model
    float min[3];
    for(int i = 0; i < 3; i++)
        min[i] = 1000000;

    float max[3];
    for(int i = 0; i < 3; i++)
        max[i] = -1000000;

    for(int k = 0; k < tc; k++)
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++) {
                min[i] = std::min(min[i], tv[k][j][i]);
                max[i] = std::max(max[i], tv[k][j][i]);
            }

    // Get the center of the model
    float center[3];
    for(int i = 0; i < 3; i++)
        center[i] = (min[i] + max[i]) / 2;

    // Get the maximum dimension of the model
    float max_dimension = 0;
    for(int i = 0; i < 3; i++)
        max_dimension = std::max(max_dimension, (max[i] - min[i]));

    // The camera looks down -Z towards the model
    // FOV is always 45 degrees, camera is moved so that the model fits.
    float camera[3];
    camera[0] = center[0];
    camera[1] = center[1];
    camera[2] = max[2] + max_dimension / .9;

    // Rather than raytrace from the camera, move the model to simplify
    // the math
    for(int k = 0; k < tc; k++) {
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                tv[k][j][i] -= camera[i];
    }

    for(int k = 0; k < tc; k++) {
        prepare_triangle(k);
    }

    printf("P2 %d %d 255\n", W, H);

    float pdx = .5 / (W / 2.0);
    float pdy = .5 / (H / 2.0);

    for(int py = 0; py < H; py++)
        for(int px = 0; px < W; px++) {

            // Make the ray
            float ray[3];
            
            ray[0] = -.5 + (px + .5) * pdx;
            ray[1] = - (-.5 + (py + .5) * pdy);
            ray[2] = -1;

            // find the closest triangle
            int triangle = -1;
            float t = F;

            for(int k = 0; k < tc; k++) {
                if(ray[0] < tb[k][0])
                    continue;
                if(ray[0] > tb[k][1])
                    continue;
                if(ray[1] < tb[k][2])
                    continue;
                if(ray[1] > tb[k][3])
                    continue;
                float t2 = intersect_triangle(k, ray);
                if(t2 < t) {
                    t = t2;
                    triangle = k;
                }
            }

            // shade the intersection
            float color;

            if(triangle == -1) {

                color = .2;

            } else {

                float facing = ray[0] * tp[triangle][0] + ray[1] * tp[triangle][1] + ray[2] * tp[triangle][2];

                float Ling = L[0] * tp[triangle][0] + L[1] * tp[triangle][1] + L[2] * tp[triangle][2];

                if(facing > 0)
                    Ling = -Ling;

                color = std::max(0.1f, Ling);
            }

            printf("%d ", (int)(color * 255));
        }
}

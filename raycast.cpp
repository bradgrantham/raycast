#include <algorithm>
#include <cmath>

// Okay, it's ugly, I get it.
//
// I didn't use classes or a vector library.  I even use loops to
// 3 everywhere.  My intent was to try to build a minimal triangle
// raycaster that I could convert simply to BASIC or FORTH or other
// simple languages without having to downconvert C++ concepts.
//
// I'm already thinking function calls is problematic, as well as
// the arrays with 3 dimensions.  I might have to inline the functions
// and use goto instead of early returns.
//
// So it's probably still not ugly enough.

const float infinity = 1000000;

const int width = 256;
const int height = 256;

const int max_triangles = 400000;

float triangle_vertices[max_triangles][3][3];
float triangle_planes[max_triangles][4];
float triangle_edges[max_triangles][3][4];

int triangle_count = 0;

float light[3] = {.577, .577, .577};

void prepare_triangle(int k)
{
    // Get edge vectors
    float e[3][3];
    for(int i = 0; i < 3; i++)
        e[0][i] = triangle_vertices[k][1][i] - triangle_vertices[k][0][i];
    for(int i = 0; i < 3; i++)
        e[1][i] = triangle_vertices[k][2][i] - triangle_vertices[k][1][i];
    for(int i = 0; i < 3; i++)
        e[2][i] = triangle_vertices[k][0][i] - triangle_vertices[k][2][i];

    // Cross product e0 and e1 to get normal
    triangle_planes[k][0] = e[0][1] * e[1][2] - e[0][2] * e[1][1];
    triangle_planes[k][1] = e[0][2] * e[1][0] - e[0][0] * e[1][2];
    triangle_planes[k][2] = e[0][0] * e[1][1] - e[0][1] * e[1][0];

    // Normalize normal
    float d = sqrt(triangle_planes[k][0] * triangle_planes[k][0] + triangle_planes[k][1] * triangle_planes[k][1] + triangle_planes[k][2] * triangle_planes[k][2]);
    for(int i = 0; i < 3; i++)
        triangle_planes[k][i] /= d;

    // Use dot product of normal and v0 to get fourth plane equation element
    triangle_planes[k][3] = triangle_planes[k][0] * triangle_vertices[k][0][0] + 
        triangle_planes[k][1] * triangle_vertices[k][0][1] + 
        triangle_planes[k][2] * triangle_vertices[k][0][2];

    // Get half-spaces representing each edge by making vector
    // perpendicular to edge and dot that with opposite vertex to get
    // plane equation
    for(int j = 0; j < 3; j++) {
        triangle_edges[k][j][0] = triangle_planes[k][1] * e[j][2] - triangle_planes[k][2] * e[j][1];
        triangle_edges[k][j][1] = triangle_planes[k][2] * e[j][0] - triangle_planes[k][0] * e[j][2];
        triangle_edges[k][j][2] = triangle_planes[k][0] * e[j][1] - triangle_planes[k][1] * e[j][0];
        triangle_edges[k][j][3] =
            triangle_edges[k][j][0] * triangle_vertices[k][j][0] +
            triangle_edges[k][j][1] * triangle_vertices[k][j][1] +
            triangle_edges[k][j][2] * triangle_vertices[k][j][2];
    }
}

float intersect_triangle(int k, float ray[3])
{
    int i;

    // find distance from ray origin (0,0,0) to plane in units of plane normal
    float factor = ray[0] * triangle_planes[k][0] + 
        ray[1] * triangle_planes[k][1] + 
        ray[2] * triangle_planes[k][2];

    // if coincident, fail
    if(factor == 0.0)
       return infinity;

    // find distance to plane in world units
    float distance = triangle_planes[k][3] / factor;

    // if intersection is behind the ray, fail
    if(distance < 0)
        return infinity;

    // calculate the point in the plane
    float point[3];
    for(int i = 0; i < 3; i++)
        point[i] = ray[i] * distance;

    for(i = 0; i < 3; i++) {
        // calculate the distance from each edge i to the point
        float edge_dot =
            triangle_edges[k][i][0] * point[0] + 
            triangle_edges[k][i][1] * point[1] + 
            triangle_edges[k][i][2] * point[2];

        // if the distance is greater than the distance to the opposite vertex, fail
        if(edge_dot < triangle_edges[k][i][3])
	   return infinity;
    }

    // succeed by returning the distance
    return distance;
}

void reflect(float incident[3], float n[3], float reflection[3])
{
    float projection = incident[0] * n[0] + incident[1] * n[1] + incident[2] * n[2];

    for(int i = 0; i < 3; i++)
        reflection[i] = incident[i] - 2.0f * projection * n[i];
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
                triangle_vertices[triangle_count][j][i] = v[j][i];
            }

        triangle_count++;
    }

    // Get the bounding box of the model
    float min[3];
    for(int i = 0; i < 3; i++)
        min[i] = 1000000;

    float max[3];
    for(int i = 0; i < 3; i++)
        max[i] = -1000000;

    for(int k = 0; k < triangle_count; k++)
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++) {
                min[i] = std::min(min[i], triangle_vertices[k][j][i]);
                max[i] = std::max(max[i], triangle_vertices[k][j][i]);
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
    for(int k = 0; k < triangle_count; k++) {
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                triangle_vertices[k][j][i] -= camera[i];
    }

    for(int k = 0; k < triangle_count; k++) {
        prepare_triangle(k);
    }

    printf("P2 %d %d 255\n", width, height);

    float pdx = .5 / (width / 2.0);
    float pdy = .5 / (height / 2.0);

    for(int py = 0; py < height; py++)
        for(int px = 0; px < width; px++) {

            // Make the ray
            float ray[3];
            
            ray[0] = -.5 + (px + .5) * pdx;
            ray[1] = - (-.5 + (py + .5) * pdy);
            ray[2] = -1;

            float d = sqrt(ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2]);

            for(int i = 0; i < 3; i++)
                ray[i] /= d;

            // find the closest triangle
            int triangle = -1;
            float t = infinity;

            for(int k = 0; k < triangle_count; k++) {
                float t2 = intersect_triangle(k, ray);
                if(t2 < t) {
                    t = t2;
                    triangle = k;
                }
            }

            // shade the intersection
            float color;

            if(triangle == -1)
                color = .2;
            else {
                float reflection[3];

                reflect(ray, triangle_planes[triangle], reflection);

                float lighting =
                    light[0] * reflection[0] + light[1] * reflection[1] + light[2] * reflection[2];

                color = std::max(0.1f, lighting);
            }

            printf("%d ", (int)(color * 255));
        }
}

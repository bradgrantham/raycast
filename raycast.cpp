#include <algorithm>
#include <cmath>

const float background = 1000000;

const int width = 256;
const int height = 256;

const int max_triangles = 65536;

float triangle_vertices[3][3][max_triangles];
float triangle_planes[4][max_triangles];
float triangle_edges[4][3][max_triangles];

int triangle_count = 0;

inline void prepare_triangle(int k)
{
    float e[3][3];
    for(int i = 0; i < 3; i++)
        e[0][i] = triangle_vertices[k][1][i] - triangle_vertices[k][0][i];
    for(int i = 0; i < 3; i++)
        e[1][i] = triangle_vertices[k][2][i] - triangle_vertices[k][1][i];
    for(int i = 0; i < 3; i++)
        e[2][i] = triangle_vertices[k][0][i] - triangle_vertices[k][2][i];

    triangle_planes[k][0] = e[0][1] * e[1][2] - e[0][2] * e[1][1];
    triangle_planes[k][1] = e[0][2] * e[1][0] - e[0][0] * e[1][2];
    triangle_planes[k][2] = e[0][0] * e[1][1] - e[0][1] * e[1][0];

    float d = sqrt(triangle_planes[k][0] * triangle_planes[k][0] + triangle_planes[k][1] * triangle_planes[k][1] + triangle_planes[k][2] * triangle_planes[k][2]);
    for(int i = 0; i < 3; i++)
        triangle_planes[k][i] /= d;

    triangle_planes[k][3] = triangle_planes[k][0] * triangle_vertices[k][0][0] + 
        triangle_planes[k][1] * triangle_vertices[k][0][1] + 
        triangle_planes[k][2] * triangle_vertices[k][0][2];

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

    float factor = ray[0] * triangle_planes[k][0] + 
        ray[1] * triangle_planes[k][1] + 
        ray[2] * triangle_planes[k][2];

    if(factor == 0.0)
       return background;

    float distance = triangle_planes[k][3] / factor;

    if(distance < 0)
        return background;

    float point[3];
    for(int i = 0; i < 3; i++)
        point[i] = ray[i] * distance;

    for(i = 0; i < 3; i++) {
        float edge_dot =
            triangle_edges[k][i][0] * point[0] + 
            triangle_edges[k][i][1] * point[1] + 
            triangle_edges[k][i][2] * point[2];

        if(edge_dot < triangle_edges[k][i][3])
	   return background;
    }

    return distance;
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

        fprintf(stderr, "%f %f %f %f %f %f %f %f %f\n",
            triangle_vertices[triangle_count][0][0], triangle_vertices[triangle_count][0][1], triangle_vertices[triangle_count][0][2],
            triangle_vertices[triangle_count][1][0], triangle_vertices[triangle_count][1][1], triangle_vertices[triangle_count][1][2],
            triangle_vertices[triangle_count][2][0], triangle_vertices[triangle_count][2][1], triangle_vertices[triangle_count][2][2]);

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
    // FOV is always 90 degrees, camera is moved so that the model fits.
    float camera[3];
    camera[0] = center[0];
    camera[1] = center[1];
    camera[2] = max[2] + max_dimension / .7;
    fprintf(stderr, "camera is at %f %f %f\n", camera[0], camera[1], camera[2]);

    // Rather than raytrace from the camera, move the model to simplify
    // the math
    for(int k = 0; k < triangle_count; k++) {
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                triangle_vertices[k][j][i] -= camera[i];
    }

    for(int k = 0; k < triangle_count; k++) {
        prepare_triangle(k);
        fprintf(stderr, "%f %f %f %f %f %f %f %f %f\n",
            triangle_vertices[k][0][0], triangle_vertices[k][0][1], triangle_vertices[k][0][2],
            triangle_vertices[k][1][0], triangle_vertices[k][1][1], triangle_vertices[k][1][2],
            triangle_vertices[k][2][0], triangle_vertices[k][2][1], triangle_vertices[k][2][2]);
        fprintf(stderr, "    plane %f %f %f %f\n",
            triangle_planes[k][0], triangle_planes[k][1], triangle_planes[k][2], triangle_planes[k][3]);
        fprintf(stderr, "    edge %f %f %f %f\n",
            triangle_edges[k][0][0], triangle_edges[k][0][1], triangle_edges[k][0][2], triangle_edges[k][0][3]);
        fprintf(stderr, "    edge %f %f %f %f\n",
            triangle_edges[k][1][0], triangle_edges[k][1][1], triangle_edges[k][1][2], triangle_edges[k][1][3]);
        fprintf(stderr, "    edge %f %f %f %f\n",
            triangle_edges[k][2][0], triangle_edges[k][2][1], triangle_edges[k][2][2], triangle_edges[k][2][3]);
    }

    printf("P2 %d %d 255\n", width, height);

    float pdx = 1.0 / (width / 2.0);
    float pdy = 1.0 / (height / 2.0);

    for(int py = 0; py < height; py++)
        for(int px = 0; px < width; px++) {

            // Make the ray
            float ray[3];
            
            ray[0] = -1 + (px + .5) * pdx;
            ray[1] = - (-1 + (py + .5) * pdy);
            ray[2] = -1;

            float d = sqrt(ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2]);

            for(int i = 0; i < 3; i++)
                ray[i] /= d;

            // find the closest triangle
            int triangle = -1;
            float t = 10000000;

            for(int k = 0; k < triangle_count; k++) {
                float t2 = intersect_triangle(k, ray);
                if(t2 < t) {
                    t = t2;
                    triangle = k;
                }
            }

            float color;

            if(triangle == -1)
                color = 0;
            else {
                color = fabs(triangle_planes[triangle][1]);
                color = 1;
            }

            printf("%d ", (int)(color * 255));
        }
}

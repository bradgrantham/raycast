#include <algorithm>
#include <cmath>

// Okay, it's ugly, I get it.
//
// I didn't use classes or a vector library.  I use loops to
// 3 everywhere for vector operations, and I even write out dot
// products.  My intent was to try to build a minimal sphere
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

const int M = 100;

float sc[M][3]; // centers
float sr[M]; // radii
float sb[M][4]; // xmin, xmax, ymin, ymax

int count = 0;

float L[3] = {.577, .577, .577};

void prepare(int k)
{
    /* 
    sb[k][0] = 2;
    sb[k][1] = -2;
    sb[k][2] = 2;
    sb[k][3] = -2;
    for(int j = 0; j < 3; j++) {
        sb[k][0] = std::min(tb[k][0], tv[k][j][0] / -tv[k][j][2]);
        sb[k][1] = std::max(tb[k][1], tv[k][j][0] / -tv[k][j][2]);
        sb[k][2] = std::min(tb[k][2], tv[k][j][1] / -tv[k][j][2]);
        sb[k][3] = std::max(tb[k][3], tv[k][j][1] / -tv[k][j][2]);
    } */
}

float intersect(int k, float r[3], float n[3])
{
    float a = r[0] * r[0] + 
        r[1] * r[1] + 
        r[2] * r[2];

    float b = 2 * (
        r[0] * (- sc[k][0]) +
        r[1] * (- sc[k][1]) +
        r[2] * (- sc[k][2])
	);

    // XXX cache this
    float c = - sr[k] * sr[k] + 
        (- sc[k][0]) * (- sc[k][0]) +
        (- sc[k][1]) * (- sc[k][1]) +
        (- sc[k][2]) * (- sc[k][2])
	;

    float det = b * b - 4 * a * c;

    if(det < 0)
        return F;

    float root = sqrt(det);

    float t1 = (- b - root) / (2 * a);
    float t2 = (- b + root) / (2 * a);

    /* "t1" is the closer intersection and "t2" is the farther intersection. */

    if(t2 < 0)
        return F;

    float t;

    if(t1 < 0)
        t = t2;
    else
        t = t1;

    float p[3];
    for(int i = 0; i < 3; i++)
        p[i] = r[i] * t;

    for(int i = 0; i < 3; i++)
        n[i] = p[i] - sc[k][i];

    // ray::normalize
    // or at least ray::length

    float d = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for(int i = 0; i < 3; i++)
        n[i] /= d;

    return t;
}

int main()
{
    float x, y, z, r;

    // Read spheres.
    while(scanf("%f %f %f %f",
        &x, &y, &z, &r) == 4) {

        sc[count][0] = x;
        sc[count][1] = y;
        sc[count][2] = z;
        sr[count] = r;
        prepare(count);

        count++;
    }

    printf("P2 %d %d 255\n", W, H);

    float pdx = .5 / (W / 2.0);
    float pdy = .5 / (H / 2.0);

    float n[3];

    for(int py = 0; py < H; py++)
        for(int px = 0; px < W; px++) {

            // Make the ray
            float r[3];
            
            r[0] = -.5 + (px + .5) * pdx;
            r[1] = - (-.5 + (py + .5) * pdy);
            r[2] = -1;

            // find the closest triangle
            int sphere = -1;
            float t = F;

            for(int k = 0; k < count; k++) {

                // if(r[0] < sb[k][0])
                    // continue;
                // if(r[0] > sb[k][1])
                    // continue;
                // if(r[1] < sb[k][2])
                    // continue;
                // if(r[1] > sb[k][3])
                    // continue;

                float n2[3];
                float t2 = intersect(k, r, n2);
                if(t2 < t) {
                    t = t2;
                    sphere = k;
                    for(int i = 0; i < 3; i++)
                        n[i] = n2[i];
                }
            }

            // shade the intersection
            float shade;

            if(sphere == -1) {

                shade = .2;

            } else {

                float facing = r[0] * n[0] + r[1] * n[1] + r[2] * n[2];

                float lighting = L[0] * n[0] + L[1] * n[1] + L[2] * n[2];

                if(facing > 0)
                    lighting = -lighting;

                shade = std::max(0.1f, lighting);
            }

            printf("%d ", (int)(shade * 255));
        }
}

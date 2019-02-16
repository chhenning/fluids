#pragma once

#include <vector>

using namespace std;
/////////
// https://github.com/CodingTrain/website/blob/master/CodingChallenges/CC_132_FluidSimulation/Processing/CC_132_FluidSimulation/Fluid.pde
// https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

struct FluidCube
{
    FluidCube(int size, float diffusion, float viscosity, float dt)
    {
        this->size = N;
        this->dt = dt;
        diff = diffusion;
        visc = viscosity;

        s = vector<float>(N * N);
        density = vector<float>(N * N);

        Vx = vector<float>(N * N);
        Vy = vector<float>(N * N);

        Vx0 = vector<float>(N * N);
        Vy0 = vector<float>(N * N);
    }

    void AddDensity(int x, int y, float amount)
    {
        int N = size;
        density[IX(x, y)] += amount;
    }

    void AddVelocity(int x, int y, float amountX, float amountY)
    {
        int N = size;
        auto index = IX(x, y);

        Vx[index] += amountX;
        Vy[index] += amountY;
    }

    void Step()
    {
        // diffuse viscosity
        diffuse(1, Vx0, Vx, visc, dt, 4, N);
        diffuse(2, Vy0, Vy, visc, dt, 4, N);

        // clean up - make sure the same amount of fluid is everywhere
        project(Vx0, Vy0, Vx, Vy, 4, N);


        advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
        advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

        // clean up
        project(Vx, Vy, Vx0, Vy0, 4, N);

        diffuse(0, s, density, diff, dt, 4, N);
        
        advect(0, density, s, Vx, Vy, dt, N);
    }

    void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter, int N)
    {
        float cRecip = 1.f / c;
        for (int k = 0; k < iter; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a * (x[IX(i + 1, j)]
                                + x[IX(i - 1, j)]
                                + x[IX(i, j + 1)]
                                + x[IX(i, j - 1)]
                                )) * cRecip;
                }
            }

            set_bnd(b, x, N);
        }
    }

    void set_bnd(int b, vector<float>& x, int N)
    {
        for (int i = 1; i < N - 1; i++) {
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
        }

        for (int j = 1; j < N - 1; j++) {
            x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
            x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
        }

        x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
        x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
        x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
    }

    void diffuse(int b, vector<float>& x, vector<float>& x0, float diff, float dt, int iter, int N)
    {
        float a = dt * diff * (N - 2) * (N - 2);
        lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
    }

    void project(vector<float>& velocX, vector<float>& velocY, vector<float>& p, vector<float>& div, int iter, int N)
    {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                    velocX[IX(i + 1, j)]
                    - velocX[IX(i - 1, j)]
                    + velocY[IX(i, j + 1)]
                    - velocY[IX(i, j - 1)]
                    ) / N;
                p[IX(i, j)] = 0;
            }
        }

        set_bnd(0, div, N);
        set_bnd(0, p, N);
        lin_solve(0, p, div, 1, 6, iter, N);

        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
                    - p[IX(i - 1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
                    - p[IX(i, j - 1)]) * N;
            }
        }

        set_bnd(1, velocX, N);
        set_bnd(2, velocY, N);
    }

    void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY, float dt, int N)
    {
        float i0, i1, j0, j1;

        float dtx = dt * (N - 2);
        float dty = dt * (N - 2);
        float dtz = dt * (N - 2);

        float s0, s1, t0, t1;
        float tmp1, tmp2, x, y;

        float Nfloat = (float)N;
        float ifloat, jfloat;
        int i, j;

        for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
            for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x = ifloat - tmp1;
                y = jfloat - tmp2;

                if (x < 0.5f) x = 0.5f;
                if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
                i0 = floorf(x);
                i1 = i0 + 1.0f;
                if (y < 0.5f) y = 0.5f;
                if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
                j0 = floorf(y);
                j1 = j0 + 1.0f;

                s1 = x - i0;
                s0 = 1.0f - s1;
                t1 = y - j0;
                t0 = 1.0f - t1;

                int i0i = (int) i0;
                int i1i = (int) i1;
                int j0i = (int) j0;
                int j1i = (int) j1;

                d[IX(i, j)] =
                    s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                    s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
            }
        }

        set_bnd(b, d, N);
    }

    int size;
    float dt;
    float diff;
    float visc;

    vector<float> s;
    vector<float> density;

    vector<float> Vx;
    vector<float> Vy;

    // previous
    vector<float> Vx0;
    vector<float> Vy0;
};


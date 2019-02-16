// Fluids.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <vector>

#define SDL_MAIN_HANDLED
#include <SDL.h>

const int N = 256;
const int Iter = 4;

//inline size_t IX(size_t x, size_t y, size_t z)
//{
//    auto ret = (N * N * z) + (N * y) + x;
//    return ret;
//}

#define IX(x, y) ((x) + (y) * N)

#include "fluids.hpp"

using namespace std;

void single_step()
{
    //const Uint64 start = SDL_GetPerformanceCounter();

    //auto cube = FluidCube(N, 0, 0, 0.1f);
    //cube.Step();

    //const Uint64 end = SDL_GetPerformanceCounter();
    //const static Uint64 freq = SDL_GetPerformanceFrequency();
    //const double seconds = (end - start) / static_cast< double >(freq);
    //cout << "Frame time: " << seconds * 1000.0 << "ms" << endl;

    //return 0;

}


int main()
{    
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        cout << "SDL_Init Error: " << SDL_GetError() << endl;
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow("Hello World!", 100, 100, 640, 480, SDL_WINDOW_SHOWN);
    if (win == nullptr) {
        cout << "SDL_CreateWindow Error: " << SDL_GetError() << endl;
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    if (ren == nullptr) {
        SDL_DestroyWindow(win);
        cout << "SDL_CreateRenderer Error: " << SDL_GetError() << endl;
        SDL_Quit();
        return 1;
    }
    
    SDL_RendererInfo info;
    SDL_GetRendererInfo(ren, &info);
    cout << "Renderer name: " << info.name << endl;
    cout << "Texture formats: " << endl;
    for (Uint32 i = 0; i < info.num_texture_formats; i++)
    {
        cout << SDL_GetPixelFormatName(info.texture_formats[i]) << endl;
    }


    const unsigned int texWidth = N;
    const unsigned int texHeight = N;
    SDL_Texture* texture = SDL_CreateTexture
    (
        ren,
        SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING,
        texWidth, texHeight
    );

    vector<unsigned char> pixels(N * N * 4, 0);

    auto cube = FluidCube(N, 0, 0, 0.1f);

    SDL_Event event;
    bool running = true;
    while (running)
    {
        const Uint64 start = SDL_GetPerformanceCounter();

        SDL_SetRenderDrawColor(ren, 0, 0, 0, SDL_ALPHA_OPAQUE);
        SDL_RenderClear(ren);

        while (SDL_PollEvent(&event))
        {
            if ((SDL_QUIT == event.type) ||
                (SDL_KEYDOWN == event.type && SDL_SCANCODE_ESCAPE == event.key.keysym.scancode))
            {
                running = false;
                break;
            }

            if (SDL_MOUSEMOTION == event.type)
            {
                SDL_MouseMotionEvent e = event.motion;
                auto mouse_x = e.x;
                auto mouse_y = e.y;

                cube.AddDensity(mouse_x, mouse_y, 100);

                printf("%d %d\n", mouse_x, mouse_y);
                //cout << mouse_x << " " << mouse_y << endl;
            }
        }

        //
        cube.Step();



        SDL_UpdateTexture
        (
            texture,
            NULL,
            &pixels[0],
            texWidth * 4
        );

        SDL_RenderCopy(ren, texture, NULL, NULL);
        SDL_RenderPresent(ren);

        const Uint64 end = SDL_GetPerformanceCounter();
        const static Uint64 freq = SDL_GetPerformanceFrequency();
        const double seconds = (end - start) / static_cast< double >(freq);
        //cout << "Frame time: " << seconds * 1000.0 << "ms" << endl;
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();



    return 0;
}


gif-h
=====

This one-header library offers a simple, very limited way to create animated GIFs directly in code.
Those looking for particular cleverness are likely to be disappointed; it's pretty much a straight-ahead
implementation of the GIF format with optional Floyd-Steinberg dithering. (It does at least use delta
encoding - only the changed portions of each frame are saved.) 

So resulting files are often quite large. The hope is that it will be handy nonetheless as a quick and easily-integrated way for programs to spit out animations.

Only RGBA8 is currently supported as an input format. (The alpha is ignored.) 

Email me : ctangora -at- gmail -dot- com

Usage:
-------------------
Create a GifWriter struct. 

Pass the struct to GifBegin() to initialize values and write the file header.

Pass frames of the animation to GifWriteFrame().

Finally, call GifEnd() to close the file handle and free memory.

    #include <gif.h>
    int main(void)
    {
        int width = 100;
        int height = 200;
        vector<uint8_t> vi1(width * height * 4, 0);   // 4 channels, RGBA
        vector<uint8_t> vi2(width * height * 4, 255);

        auto fileName = "bwgif.gif";
        int delay = 100; // 100ms
        GifWriter g;
        GifBegin(&g, fileName, width, height, delay);
        GifWriteFrame(&g, vi1.data(), width, height, delay);
        GifWriteFrame(&g, vi2.data(), width, height, delay);
        GifEnd(&g);

        return 0;
    }

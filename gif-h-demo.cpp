//
// gif-h-demo.cpp
// by Charlie Tangora
// Public domain.
// Email me : ctangora -at- gmail -dot- com
//
// Shows an example usage of gif.h
//

#include "gif.h"

#include <math.h>

const int width = 256;
const int height = 256;
uint8_t image[ width * height * 4 ];

void SetPixel( int xx, int yy, uint8_t red, uint8_t grn, uint8_t blu )
{
    uint8_t* pixel = &image[(yy*width+xx)*4];
    pixel[0] = red;
    pixel[1] = blu;
    pixel[2] = grn;
    pixel[3] = 255;  // no alpha for this demo
}

void SetPixelFloat( int xx, int yy, float fred, float fgrn, float fblu )
{
    // convert float to unorm
    uint8_t red = (uint8_t)roundf( 255.0f * fred );
    uint8_t grn = (uint8_t)roundf( 255.0f * fgrn );
    uint8_t blu = (uint8_t)roundf( 255.0f * fblu );
    
    SetPixel( xx, yy, red, grn, blu );
}

int main(int argc, const char * argv[])
{
    const char* filename = "./MyGif.gif";
    if( argc > 1 )
    {
        filename = argv[1];
    }
    
    // Create a gif
    GifWriter writer = {};
    GifBegin( &writer, filename, width, height, 2, 8, true );
    
    for( int frame=0; frame<256; ++frame )
    {
        
        // Make an image, somehow
        // this is the default shadertoy - credit to shadertoy.com
        float tt = frame * 3.14159f * 2 / 255.0f;
        for( int yy=0; yy<height; ++yy )
        {
            for( int xx=0; xx<width; ++xx )
            {
                float fx = xx / (float)width;
                float fy = yy / (float)height;
                
                float red = 0.5f + 0.5f * cosf(tt+fx);
                float grn = 0.5f + 0.5f * cosf(tt+fy+2.f);
                float blu = 0.5f + 0.5f * cosf(tt+fx+4.f);
                
                SetPixelFloat( xx, yy, red, grn, blu );
            }
        }
        
        
        // Write the frame to the gif
        printf( "Writing frame %d...\n", frame );
        GifWriteFrame( &writer, image, width, height, 2, 8, true );
    }
    
    // Write EOF
    GifEnd( &writer );
    
    return 0;
}

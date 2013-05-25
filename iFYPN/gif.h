//
// gif.h
// by Charlie Tangora
//
// This file offers a simple, very limited way to create animated GIFs directly in code.
// There's not much effort to optimize colors to improve compression, and only a 256-color
// palette is supported, so resulting files are quite large. But it's still handy to be able
// to spit out GIFs on command.
//
// Input images must be RGBA8 formatted.
//
// USAGE:
// Create a GifWriter struct. Pass it to GifBegin() with the first frame.
// Pass subsequent frames to GifContinue().
// Finally, call GifEnd() to close the file handle and free memory.
//

#ifndef gif_h
#define gif_h

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"

const int kGifTransIndex = 0;

struct GifPalette
{
    uint8_t r[256];
    uint8_t g[256];
    uint8_t b[256];
    
    uint8_t treeSplitElt[255];
    uint8_t treeSplit[255];
};

// walks the k-d tree to pick the palette entry for a desired color.
// Takes as in/out parameters the current best color and its error -
// only changes them if it finds a better color in its subtree.
// this is the major hotspot in the code at the moment.
void GifGetClosestPaletteColor(GifPalette* pPal, int r, int g, int b, int& bestInd, int& bestDiff, int treeRoot = 1)
{
    // base case, reached the bottom of the tree
    if(treeRoot > 255)
    {
        int ind = treeRoot-256;
        if(ind == kGifTransIndex) return;
        
        // check whether this color is better than the current winner
        int r_err = r - ((int32_t)pPal->r[ind]);
        int g_err = g - ((int32_t)pPal->g[ind]);
        int b_err = b - ((int32_t)pPal->b[ind]);
        int diff = abs(r_err)+abs(g_err)+abs(b_err);
        
        if(diff < bestDiff)
        {
            bestInd = ind;
            bestDiff = diff;
        }
        
        return;
    }
    
    // take the appropriate color (r, g, or b) for this node of the k-d tree
    int comps[3]; comps[0] = r; comps[1] = g; comps[2] = b;
    int splitComp = comps[pPal->treeSplitElt[treeRoot]];
    
    uint32_t splitPos = pPal->treeSplit[treeRoot];
    if(splitPos > splitComp)
    {
        // check the left subtree
        GifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2);
        if( bestDiff > splitPos - splitComp )
        {
            // there may be a better value in the right subtree, check that too
            GifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2+1);
        }
    }
    else
    {
        GifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2+1);
        if( bestDiff > splitComp - splitPos )
        {
            GifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2);
        }
    }
}

void GifSwapPixels(uint8_t* image, int pixA, int pixB)
{
    uint8_t rA = image[pixA*4];
    uint8_t gA = image[pixA*4+1];
    uint8_t bA = image[pixA*4+2];
    uint8_t aA = image[pixA*4+3];
    
    uint8_t rB = image[pixB*4];
    uint8_t gB = image[pixB*4+1];
    uint8_t bB = image[pixB*4+2];
    uint8_t aB = image[pixA*4+3];
    
    image[pixA*4] = rB;
    image[pixA*4+1] = gB;
    image[pixA*4+2] = bB;
    image[pixA*4+3] = aB;
    
    image[pixB*4] = rA;
    image[pixB*4+1] = gA;
    image[pixB*4+2] = bA;
    image[pixB*4+3] = aA;
}

// just the partition operation from quicksort
int GifPartition(uint8_t* image, const int left, const int right, const int elt, int pivotIndex)
{
    const int pivotValue = image[(pivotIndex)*4+elt];
    GifSwapPixels(image, pivotIndex, right-1);
    int storeIndex = left;
    bool split = 0;
    for(int ii=left; ii<right-1; ++ii)
    {
        int arrayVal = image[ii*4+elt];
        if( arrayVal < pivotValue )
        {
            GifSwapPixels(image, ii, storeIndex);
            ++storeIndex;
        }
        else if( arrayVal == pivotValue )
        {
            if(split)
            {
                GifSwapPixels(image, ii, storeIndex);
                ++storeIndex;
            }
            split = !split;
        }
    }
    GifSwapPixels(image, storeIndex, right-1);
    return storeIndex;
}

// Perform an incomplete sort, finding all elements above and below the desired median
void GifPartitionByMedian(uint8_t* image, int left, int right, int com, int neededCenter)
{
    if(left < right-1)
    {
        int pivotIndex = left + (right-left)/2;
    
        pivotIndex = GifPartition(image, left, right, com, pivotIndex);
        
        // Only "sort" the section of the array that contains the median
        if(pivotIndex > neededCenter)
            GifPartitionByMedian(image, left, pivotIndex, com, neededCenter);
        
        if(pivotIndex < neededCenter)
            GifPartitionByMedian(image, pivotIndex+1, right, com, neededCenter);
    }
}

// Builds a palette by creating a balanced k-d tree of all pixels in the image
void GifSplitPalette(uint8_t* image, int numPixels, int firstElt, int lastElt, int splitElt, int splitDist, int treeNode, GifPalette* pal)
{
    if(lastElt <= firstElt)
        return;
    
    // base case, bottom of the tree
    if(lastElt == firstElt+1)
    {
        // take the average of all colors in this subcube
        // note: for color index 1 or 255 (the ends)
        // it might be better to use the max/min color
        // in the subcube - otherwise Floyd-Steinberg
        // dithering can get confused.
        // Doesn't seem to be a big issue so far, though.
        
        int r=0, g=0, b=0;
        for(int ii=0; ii<numPixels; ++ii)
        {
            r += image[ii*4+0];
            g += image[ii*4+1];
            b += image[ii*4+2];
        }
        
        r += numPixels / 2;  // round to nearest
        g += numPixels / 2;
        b += numPixels / 2;
        
        r /= numPixels;
        g /= numPixels;
        b /= numPixels;
        
        pal->r[firstElt] = r;
        pal->g[firstElt] = g;
        pal->b[firstElt] = b;
        
        return;
    }
    
    // Find the axis with the largest range
    int minR = 255, maxR = 0;
    int minG = 255, maxG = 0;
    int minB = 255, maxB = 0;
    for(int ii=0; ii<numPixels; ++ii)
    {
        int r = image[ii*4+0];
        int g = image[ii*4+1];
        int b = image[ii*4+2];
        
        if(r > maxR) maxR = r;
        if(r < minR) minR = r;
        
        if(g > maxG) maxG = g;
        if(g < minG) minG = g;
        
        if(b > maxB) maxB = b;
        if(b < minB) minB = b;
    }
    
    int rRange = maxR - minR;
    int gRange = maxG - minG;
    int bRange = maxB - minB;
    
    // and split along that axis. (incidentally, this means this isn't a "proper" k-d tree but I don't know what else to call it)
    int splitCom = 1;
    if(bRange > gRange) splitCom = 2;
    if(rRange > bRange && rRange > gRange) splitCom = 0;
    
    int subPixelsA = numPixels * (splitElt - firstElt) / (lastElt - firstElt);
    int subPixelsB = numPixels-subPixelsA;
    
    GifPartitionByMedian(image, 0, numPixels, splitCom, subPixelsA);
    
    pal->treeSplitElt[treeNode] = splitCom;
    pal->treeSplit[treeNode] = image[subPixelsA*4+splitCom];
    
    GifSplitPalette(image,              subPixelsA, firstElt, splitElt, splitElt-splitDist, splitDist/2, treeNode*2,   pal);
    GifSplitPalette(image+subPixelsA*4, subPixelsB, splitElt, lastElt,  splitElt+splitDist, splitDist/2, treeNode*2+1, pal);
}

// Creates a palette by placing all the image pixels in a k-d tree and then averaging the blocks at the bottom.
// This is known as the "modified median split" technique
void GifMakePalette( uint8_t* image, uint32_t width, uint32_t height, GifPalette* pPal )
{
    // SplitPalette is destructive (it sorts the pixels by color) so
    // we must create a copy of the image for it to destroy
    int imageSize = width*height*4*sizeof(uint8_t);
    uint8_t* destroyableImage = (uint8_t*)malloc(imageSize);
    memcpy(destroyableImage, image, imageSize);
    
    GifSplitPalette(destroyableImage, width*height, 1, 256, 128, 64, 1, pPal);
    
    free(destroyableImage);
    
    // add the bottom node for the transparency index
    pPal->treeSplit[128] = 0;
    pPal->treeSplitElt[128] = 0;
    
    pPal->r[0] = pPal->g[0] = pPal->b[0] = 0;
}

int GifIMax(int l, int r)
{
    return l>r?l:r;
}

// Implements Floyd-Steinberg dithering, writes palette value to alpha
void GifDitherImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, uint32_t width, uint32_t height, GifPalette* pPal )
{
    uint32_t numPixels = width*height;
    
    // quantPixels initially holds color*256 for all pixels
    // The extra 8 bits of precision allow for sub-single-color error values
    // to be propagated
    int32_t* quantPixels = (int32_t*)malloc(sizeof(int32_t)*numPixels*4);
    
    for( uint32_t ii=0; ii<numPixels*4; ++ii )
    {
        uint8_t pix = nextFrame[ii];
        uint16_t pix16 = (uint16_t)(pix) << 8;
        quantPixels[ii] = pix16;
    }
    
    for( uint32_t yy=0; yy<height; ++yy )
    {
        for( uint32_t xx=0; xx<width; ++xx )
        {
            int32_t* nextPix = quantPixels + 4*(yy*width+xx);
            const uint8_t* lastPix = lastFrame? lastFrame + 4*(yy*width+xx) : NULL;
            
            // initialize the search with the transparent value - if it's best, we write out a transparent pixel
            int32_t r_err = lastFrame? ((int32_t)nextPix[0] >> 8) - ((int32_t)lastPix[0]) : 1000000;
            int32_t g_err = lastFrame? ((int32_t)nextPix[1] >> 8) - ((int32_t)lastPix[1]) : 1000000;
            int32_t b_err = lastFrame? ((int32_t)nextPix[2] >> 8) - ((int32_t)lastPix[2]) : 1000000;
            
            int32_t bestDiff = abs(r_err)+abs(g_err)+abs(b_err);
            int32_t bestInd = kGifTransIndex;
            
            // Search the palete
            int32_t rr = ((int32_t)nextPix[0] + 128) >> 8;
            int32_t gg = ((int32_t)nextPix[1] + 128) >> 8;
            int32_t bb = ((int32_t)nextPix[2] + 128) >> 8;
            GifGetClosestPaletteColor(pPal, rr, gg, bb, bestInd, bestDiff);
            
            // Write the result to the temp buffer
            if(bestInd != kGifTransIndex)
            {
                r_err = ((int32_t)nextPix[0]) - (((int32_t)pPal->r[bestInd]) << 8);
                g_err = ((int32_t)nextPix[1]) - (((int32_t)pPal->g[bestInd]) << 8);
                b_err = ((int32_t)nextPix[2]) - (((int32_t)pPal->b[bestInd]) << 8);
                
                nextPix[0] = pPal->r[bestInd];
                nextPix[1] = pPal->g[bestInd];
                nextPix[2] = pPal->b[bestInd];
            }
            
            if( lastFrame &&
               lastPix[0] == nextPix[0] &&
               lastPix[1] == nextPix[1] &&
               lastPix[2] == nextPix[2] )
            {
                nextPix[3] = kGifTransIndex;
            }
            else
            {
                nextPix[3] = bestInd;
            }
            
            // Propagate the error to the four adjacent locations
            // that we haven't touched yet
            int quantloc_7 = (yy*width+xx+1);
            int quantloc_3 = (yy*width+width+xx-1);
            int quantloc_5 = (yy*width+width+xx);
            int quantloc_1 = (yy*width+width+xx+1);
            
            if(quantloc_7 < numPixels)
            {
                int32_t* pix7 = quantPixels+4*quantloc_7;
                pix7[0] += GifIMax( -pix7[0], r_err * 7 / 16 );
                pix7[1] += GifIMax( -pix7[1], g_err * 7 / 16 );
                pix7[2] += GifIMax( -pix7[2], b_err * 7 / 16 );
            }
            
            if(quantloc_3 < numPixels)
            {
                int32_t* pix3 = quantPixels+4*quantloc_3;
                pix3[0] += GifIMax( -pix3[0], r_err * 3 / 16 );
                pix3[1] += GifIMax( -pix3[1], g_err * 3 / 16 );
                pix3[2] += GifIMax( -pix3[2], b_err * 3 / 16 );
            }
            
            if(quantloc_5 < numPixels)
            {
                int32_t* pix5 = quantPixels+4*quantloc_5;
                pix5[0] += GifIMax( -pix5[0], r_err * 5 / 16 );
                pix5[1] += GifIMax( -pix5[1], g_err * 5 / 16 );
                pix5[2] += GifIMax( -pix5[2], b_err * 5 / 16 );
            }
            
            if(quantloc_1 < numPixels)
            {
                int32_t* pix1 = quantPixels+4*quantloc_1;
                pix1[0] += GifIMax( -pix1[0], r_err / 16 );
                pix1[1] += GifIMax( -pix1[1], g_err / 16 );
                pix1[2] += GifIMax( -pix1[2], b_err / 16 );
            }
        }
    }
    
    // Copy the palettized result to the output buffer
    for( uint32_t ii=0; ii<numPixels*4; ++ii )
    {
        outFrame[ii] = quantPixels[ii];
    }
    
    free(quantPixels);
}

// Picks palette colors for the image using simple thresholding, no dithering
void GifThresholdImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, uint32_t width, uint32_t height, GifPalette* pPal )
{
    uint32_t numPixels = width*height;
    for( uint32_t ii=0; ii<numPixels; ++ii )
    {
        int32_t bestDiff = 1000000;
        int32_t bestInd = 1;
        
        // if a previous color is available, set the current best color to "transparent"
        if(lastFrame)
        {
            int32_t r_err = ((int32_t)lastFrame[0]) - ((int32_t)nextFrame[0]);
            int32_t g_err = ((int32_t)lastFrame[1]) - ((int32_t)nextFrame[1]);
            int32_t b_err = ((int32_t)lastFrame[2]) - ((int32_t)nextFrame[2]);
            
            bestDiff = abs(r_err)+abs(g_err)+abs(b_err);
            bestInd = kGifTransIndex;
        }
        
        GifGetClosestPaletteColor(pPal, nextFrame[0], nextFrame[1], nextFrame[2], bestInd, bestDiff);
        ALWAYS_ASSERT(lastFrame || bestInd != kGifTransIndex);
        
        // Write the resulting color to the output buffer
        if(bestInd == kGifTransIndex)
        {
            outFrame[0] = lastFrame[0];
            outFrame[1] = lastFrame[1];
            outFrame[2] = lastFrame[2];
        }
        else
        {
            outFrame[0] = pPal->r[bestInd];
            outFrame[1] = pPal->g[bestInd];
            outFrame[2] = pPal->b[bestInd];
        }
        outFrame[3] = bestInd;
        
        if(lastFrame) lastFrame += 4;
        outFrame += 4;
        nextFrame += 4;
    }
}

// Simple structure to write out the LZW-compressed portion of the image
// one bit at a time
struct GifBitStatus
{
    uint8_t bitIndex;  // how many bits in the partial byte written so far
    uint8_t byte;      // current partial byte
    
    uint32_t chunkIndex;
    uint8_t chunk[256];   // bytes are written in here until we have 256 of them, then written to the file
};

// insert a single bit
void GifWriteBit( GifBitStatus& stat, uint32_t bit )
{
    bit = bit & 1;
    bit = bit << stat.bitIndex;
    stat.byte |= bit;
    
    ++stat.bitIndex;
    if( stat.bitIndex > 7 )
    {
        // move the newly-finished byte to the chunk buffer 
        stat.chunk[stat.chunkIndex++] = stat.byte;
        // and start a new byte
        stat.bitIndex = 0;
        stat.byte = 0;
    }
}

// write all bytes so far to the file
void GifWriteChunk( FILE* f, GifBitStatus& stat )
{
    fputc(stat.chunkIndex, f);
    fwrite(stat.chunk, 1, stat.chunkIndex, f);
    
    stat.bitIndex = 0;
    stat.byte = 0;
    stat.chunkIndex = 0;
}

void GifWriteCode( FILE* f, GifBitStatus& stat, uint32_t code, uint32_t length )
{
    for( uint32_t ii=0; ii<length; ++ii )
    {
        GifWriteBit(stat, code);
        code = code >> 1;
        
        if( stat.chunkIndex == 255 )
        {
            GifWriteChunk(f, stat);
        }
    }
}

// The LZW dictionary is a 256-ary tree constructed as the file is encoded,
// this is one node
struct GifLzwNode
{
    uint16_t m_next[256];
};

// write a 256-color (8-bit) image palette to the file
void GifWritePalette( const GifPalette* pPal, FILE* f )
{
    fputc(0, f);  // first color: transparency
    fputc(0, f);
    fputc(0, f);
    
    for(uint32_t ii=1; ii<256; ++ii)
    {
        uint32_t r = pPal->r[ii];
        uint32_t g = pPal->g[ii];
        uint32_t b = pPal->b[ii];
        
        fputc(r, f);
        fputc(g, f);
        fputc(b, f);
    }
}

// write the image header, LZW-compress and write out the image
void GifWriteLzwImage(FILE* f, uint8_t* image, uint32_t left, uint32_t top,  uint32_t width, uint32_t height, uint32_t delay, GifPalette* pPal)
{
    // graphics control extension
    fputc(0x21, f);
    fputc(0xf9, f);
    fputc(0x04, f);
    fputc(0x05, f); // dispose by leaving in place, has transparency
    fputc(delay & 0xff, f);
    fputc((delay >> 8) & 0xff, f);
    fputc(kGifTransIndex, f); // transparent color index
    fputc(0, f);
    
    fputc(0x2c, f); // image descriptor block
    
    fputc(left & 0xff, f);           // corner of image in canvas space
    fputc((left >> 8) & 0xff, f);
    fputc(top & 0xff, f);
    fputc((top >> 8) & 0xff, f);
    
    fputc(width & 0xff, f);          // width and height of image
    fputc((width >> 8) & 0xff, f);
    fputc(height & 0xff, f);
    fputc((height >> 8) & 0xff, f);
    
    //fputc(0, f); // no local color table, no transparency
    //fputc(0x80, f); // no local color table, but transparency
    
    fputc(0x87, f); // local color table present, 256 entries
    GifWritePalette(pPal, f);
    
    fputc(8, f); // min code size 8 bits
    
    GifLzwNode* codetree = (GifLzwNode*)malloc(sizeof(GifLzwNode)*4096);
    
    bzero(codetree, sizeof(GifLzwNode)*4096);
    int32_t curCode = -1;
    uint32_t codeSize = 9;
    uint32_t maxCode = 257;
    
    GifBitStatus stat;
    stat.byte = 0;
    stat.bitIndex = 0;
    stat.chunkIndex = 0;
    
    GifWriteCode(f, stat, 256, codeSize);  // start with a fresh LZW dictionary
    
    for(int32_t yy=0; yy<height; ++yy)
    {
        for(int32_t xx=0; xx<width; ++xx)
        {
            uint8_t nextValue = image[(yy*width+xx)*4+3];
            
            // "loser mode" - no compression, every single code is followed immediately by a clear
            //WriteCode( f, stat, nextValue, codeSize );
            //WriteCode( f, stat, 256, codeSize );
            
            if( curCode < 0 )
            {
                // first value in a new run
                curCode = nextValue;
            }
            else if( codetree[curCode].m_next[nextValue] )
            {
                // current run already in the dictionary
                curCode = codetree[curCode].m_next[nextValue];
            }
            else
            {
                // finish the current run, write a code
                GifWriteCode( f, stat, curCode, codeSize );
                
                // insert the new run into the dictionary
                codetree[curCode].m_next[nextValue] = ++maxCode;
                
                if( maxCode >= (1 << codeSize) )
                {
                    // dictionary entry count has broken a size barrier,
                    // we need more bits for codes
                    codeSize++;
                }
                if( maxCode == 4095 )
                {
                    // the dictionary is full, clear it out and begin anew
                    GifWriteCode(f, stat, 256, codeSize); // clear tree
                    
                    bzero(codetree, sizeof(GifLzwNode)*4096);
                    curCode = -1;
                    codeSize = 9;
                    maxCode = 257;
                }
                
                curCode = nextValue;
            }
        }
    }
    
    // compression footer
    GifWriteCode( f, stat, curCode, codeSize );
    GifWriteCode( f, stat, 256, codeSize );
    GifWriteCode( f, stat, 257, 9 );
    
    // write out the last partial chunk
    while( stat.bitIndex ) GifWriteBit(stat, 0);
    if( stat.chunkIndex ) GifWriteChunk(f, stat);
    
    fputc(0, f); // image block terminator
    
    free(codetree);
}

struct GifWriter
{
    FILE* f;
    uint8_t* oldImage;
};

// Creates a gif file, writes out the first frame.
// The input GIFWriter is assumed to be uninitialized.
// The delay value is the time between frames in hundredths of a second - note that not all viewers pay much attention to this value.
bool GifBegin( GifWriter* writer, const char* filename, uint8_t* image, uint32_t width, uint32_t height, uint32_t delay, bool dither = false )
{
    writer->f = fopen(filename, "wb");
    if(!writer->f) return false;
    
    // allocate 
    writer->oldImage = (uint8_t*)malloc(width*height*4);
    
    GifPalette pal;
    GifMakePalette(image, width, height, &pal);
    
    if(dither)
        GifDitherImage(NULL, image, writer->oldImage, width, height, &pal);
    else
        GifThresholdImage(NULL, image, writer->oldImage, width, height, &pal);
    
    fputs("GIF89a", writer->f);
    
    // screen descriptor
    fputc(width & 0xff, writer->f);
    fputc((width >> 8) & 0xff, writer->f);
    fputc(height & 0xff, writer->f);
    fputc((height >> 8) & 0xff, writer->f);
    
    fputc(0xf7, writer->f);  // there is an unsorted global color table of 256 entries
    fputc(0, writer->f);     // background color
    fputc(0, writer->f);     // pixels are square (we need to specify this because it's 1989)
    
    // now the global palette
    GifWritePalette(&pal, writer->f);
    
    if( delay != 0 )
    {
        // animation header
        fputc(0x21, writer->f); // extension
        fputc(0xff, writer->f); // application specific
        fputc(11, writer->f); // length 11
        fputs("NETSCAPE2.0", writer->f); // yes, really
        fputc(3, writer->f); // 3 bytes of NETSCAPE2.0 data
        
        fputc(1, writer->f); // JUST BECAUSE
        fputc(0, writer->f); // loop infinitely (byte 0)
        fputc(0, writer->f); // loop infinitely (byte 1)
        
        fputc(0, writer->f); // block terminator
    }
    
    GifWriteLzwImage(writer->f, writer->oldImage, 0, 0, width, height, delay, &pal);
    
    return true;
}

// Writes out a new frame to a GIF in progress.
// The GIFWriter should have been created by GIFBegin.
bool GifContinue( GifWriter* writer, uint8_t* image, uint32_t width, uint32_t height, uint32_t delay, bool dither = false )
{
    if(!writer->f) return false;
    
    GifPalette pal;
    GifMakePalette(image, width, height, &pal);
    
    if(dither)
        GifDitherImage(writer->oldImage, image, writer->oldImage, width, height, &pal);
    else
        GifThresholdImage(writer->oldImage, image, writer->oldImage, width, height, &pal);
    
    GifWriteLzwImage(writer->f, writer->oldImage, 0, 0, width, height, delay, &pal);
    
    return true;
}

// Writes the EOF code, closes the file handle, and frees temp memory used by a GIF.
// Many if not most viewers will still display a GIF properly if the EOF code is missing,
// but it's still a good idea to write it out.
bool GifEnd( GifWriter* writer )
{
    if(!writer->f) return false;
    
    fputc(0x3b, writer->f); // end of file
    fclose(writer->f);
    free(writer->oldImage);
    
    writer->f = NULL;
    writer->oldImage = NULL;
    
    return true;
}

#endif

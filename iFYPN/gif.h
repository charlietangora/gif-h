//
//  gif.h
//  OSXGLEssentials
//
//  Created by Charles Tangora on 8/22/12.
//
//

#ifndef gif_h
#define gif_h

#include <stdio.h>
#include <stdint.h>

struct Palette
{
    uint32_t r[256];
    uint32_t g[256];
    uint32_t b[256];
};

void SwapPixels(uint8_t* image, int pixA, int pixB)
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

int Partition(uint8_t* image, const int left, const int right, const int elt, int pivotIndex)
{
    const int pivotValue = image[(pivotIndex)*4+elt];
    SwapPixels(image, pivotIndex, right-1);
    int storeIndex = left;
    bool split = 0;
    for(int ii=left; ii<right-1; ++ii)
    {
        int arrayVal = image[ii*4+elt];
        if( arrayVal < pivotValue )
        {
            SwapPixels(image, ii, storeIndex);
            ++storeIndex;
        }
        else if( arrayVal == pivotValue )
        {
            if(split)
            {
                SwapPixels(image, ii, storeIndex);
                ++storeIndex;
            }
            split = !split;
        }
    }
    SwapPixels(image, storeIndex, right-1);
    return storeIndex;
}

// Perform an incomplete sort, finding all elements above and below the desired median
void PartitionByMedian(uint8_t* image, int left, int right, int elt, int neededCenter)
{
    if(left < right-1)
    {
        int pivotIndex = left + (right-left)/2;
    
        pivotIndex = Partition(image, left, right, elt, pivotIndex);
        
        if(pivotIndex > neededCenter)
            PartitionByMedian(image, left, pivotIndex, elt, neededCenter);
        
        if(pivotIndex < neededCenter)
            PartitionByMedian(image, pivotIndex+1, right, elt, neededCenter);
    }
}

void SplitPalette(uint8_t* image, int numPixels, int firstElt, int numElts, Palette* pal)
{
    if(numElts == 1)
    {
        int r=0, g=0, b=0;
        for(int ii=0; ii<numPixels; ++ii)
        {
            r += image[ii*4+0];
            g += image[ii*4+1];
            b += image[ii*4+2];
        }
        
        r += numPixels / 2;
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
    
    int splitElt = 1;
    if(bRange > gRange) splitElt = 2;
    if(rRange > bRange && rRange > gRange) splitElt = 0;
    
    int subEltsA = numElts/2;
    int subEltsB = numElts-subEltsA;
    
    int subPixelsA = numPixels * subEltsA / numElts;
    int subPixelsB = numPixels-subPixelsA;
    
    PartitionByMedian(image, 0, numPixels, splitElt, subPixelsA);
    
    SplitPalette(image,              subPixelsA, firstElt,          subEltsA, pal);
    SplitPalette(image+subPixelsA*4, subPixelsB, firstElt+subEltsA, subEltsB, pal);
}

void MakePalette( uint8_t* image, uint32_t width, uint32_t height, Palette* pPal )
{
    pPal->r[0] = pPal->g[0] = pPal->b[0] = 0;
    pPal->r[1] = pPal->g[1] = pPal->b[1] = 0;
    
    pPal->r[255] = pPal->g[255] = pPal->b[255] = 255;
    pPal->r[254] = pPal->g[254] = pPal->b[254] = 0;
    pPal->r[253] = pPal->g[253] = pPal->b[253] = 0;
    pPal->r[252] = pPal->g[252] = pPal->b[252] = 0;
    
    for(int ii=0; ii<width*height; ++ii)
    {
        int r = image[ii*4+0];
        int g = image[ii*4+1];
        int b = image[ii*4+2];
        
        if( r > pPal->r[254] || (r == pPal->r[254] && (g > pPal->g[254] || b > pPal->b[254])) )
        {
            pPal->r[254] = r;
            pPal->g[254] = g;
            pPal->b[254] = b;
        }
        
        if( g > pPal->g[253] || (g == pPal->g[253] && (r > pPal->r[253] || b > pPal->b[253])) )
        {
            pPal->r[253] = r;
            pPal->g[253] = g;
            pPal->b[253] = b;
        }
        
        if( b > pPal->b[252] || (b == pPal->b[252] && (r > pPal->r[252] || b > pPal->b[252])) )
        {
            pPal->r[252] = r;
            pPal->g[252] = g;
            pPal->b[252] = b;
        }
    }
    
    int imageSize = width*height*4*sizeof(uint8_t);
    uint8_t* destroyableImage = (uint8_t*)malloc(imageSize);
    memcpy(destroyableImage, image, imageSize);
    
    SplitPalette(destroyableImage, width*height, 2, 250, pPal);
    
    free(destroyableImage);
}

struct BitStatus
{
    uint8_t bitIndex;
    uint8_t byte;
    
    uint32_t chunkIndex;
    uint8_t chunk[256];
};

void WriteBit( BitStatus& stat, uint32_t bit )
{
    bit = bit & 1;
    bit = bit << stat.bitIndex;
    stat.byte |= bit;
    
    ++stat.bitIndex;
    if( stat.bitIndex > 7 )
    {
        stat.chunk[stat.chunkIndex++] = stat.byte;
        stat.bitIndex = 0;
        stat.byte = 0;
    }
}

void WriteImageChunk( FILE* f, BitStatus& stat )
{
    ASSERT(stat.bitIndex == 0);
    ASSERT(stat.chunkIndex < 256);
    
    fputc(stat.chunkIndex, f);
    fwrite(stat.chunk, 1, stat.chunkIndex, f);
    
    stat.bitIndex = 0;
    stat.byte = 0;
    stat.chunkIndex = 0;
}

void WriteCode( FILE* f, BitStatus& stat, uint32_t code, uint32_t length )
{
    for( uint32_t ii=0; ii<length; ++ii )
    {
        WriteBit(stat, code);
        code = code >> 1;
        
        if( stat.chunkIndex == 255 )
        {
            WriteImageChunk(f, stat);
        }
    }
}

struct lzwNode
{
    uint32_t m_code;
    lzwNode* m_children[256];
};

lzwNode* InitLzwTree()
{
    lzwNode* root = (lzwNode*)malloc(sizeof(lzwNode));
    for(int ii=0; ii<256; ++ii)
    {
        lzwNode* child = (lzwNode*)malloc(sizeof(lzwNode));
        bzero(child, sizeof(lzwNode));
        child->m_code = ii;
        root->m_children[ii] = child;
    }
    return root;
}

void DeleteLzwTree( lzwNode* node )
{
    if( !node ) return;
    
    for( uint32_t ii=0; ii<256; ++ii )
        DeleteLzwTree(node->m_children[ii]);
    
    free(node);
}

void WritePalette( const Palette* pPal, FILE* f )
{
    fputc(0, f);  // first and second colors: always black
    fputc(0, f);
    fputc(0, f);
    
    fputc(0, f);
    fputc(0, f);
    fputc(0, f);
    
    for(uint32_t ii=2; ii<255; ++ii)
    {
        uint32_t r = pPal->r[ii];
        uint32_t g = pPal->g[ii];
        uint32_t b = pPal->b[ii];
        
        fputc(r, f);
        fputc(g, f);
        fputc(b, f);
    }
    
    fputc(255, f);  // last color: always white
    fputc(255, f);
    fputc(255, f);
}

void WriteLzwImage(FILE* f, uint8_t* image, uint8_t* oldImage, uint32_t left, uint32_t top,  uint32_t width, uint32_t height, uint32_t delay, Palette& pal)
{
    // graphics control extension
    fputc(0x21, f);
    fputc(0xf9, f);
    fputc(0x04, f);
    fputc(0x05, f); // dispose by leaving in place, has transparency
    fputc(delay & 0xff, f);
    fputc((delay >> 8) & 0xff, f);
    fputc(0x01, f); // transparent color index 1
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
    WritePalette(&pal, f);
    
    fputc(8, f); // min code size 8 bits
    
    lzwNode* codetree = InitLzwTree();
    lzwNode* curNode = codetree;
    uint32_t codeSize = 9;
    uint32_t maxCode = 257;
    
    BitStatus stat;
    stat.byte = 0;
    stat.bitIndex = 0;
    stat.chunkIndex = 0;
    
    WriteCode(f, stat, 256, codeSize);
    
    for(int32_t yy=0; yy<height; ++yy)
    {
        for(int32_t xx=0; xx<width; ++xx)
        {
            uint8_t nextValue = image[(yy*width+xx)*4+3];
            
            // "loser mode" - no compression, every single code is followed immediately by a clear
            //WriteCode( f, stat, nextValue, codeSize );
            //WriteCode( f, stat, 256, codeSize );
            
            if( curNode->m_children[nextValue] )
            {
                curNode = curNode->m_children[nextValue];
            }
            else
            {
                WriteCode( f, stat, curNode->m_code, codeSize );
                
                lzwNode* child = (lzwNode*)malloc(sizeof(lzwNode));
                bzero(child, sizeof(lzwNode));
                maxCode++;
                child->m_code = maxCode;
                curNode->m_children[nextValue] = child;
                
                if( maxCode >= (1 << codeSize) )
                {
                    codeSize++;
                }
                if( maxCode == 4095 )
                {
                    ASSERT(codeSize == 12);
                    WriteCode(f, stat, 256, codeSize); // clear tree
                    
                    DeleteLzwTree(codetree);
                    codetree = InitLzwTree();
                    curNode = codetree;
                    codeSize = 9;
                    maxCode = 257;
                }
                
                curNode = codetree->m_children[nextValue];
            }
        }
    }
    
    WriteCode( f, stat, curNode->m_code, codeSize );
    WriteCode( f, stat, 256, codeSize );
    WriteCode( f, stat, 257, 9 );
    while( stat.bitIndex ) WriteBit(stat, 0);
    if( stat.chunkIndex ) WriteImageChunk(f, stat);
    DeleteLzwTree(codetree);
    fputc(0, f); // image block terminator
}

void SetTransparency( uint8_t* lastFrame, uint8_t* nextFrame, uint32_t width, uint32_t height )
{
    uint32_t numPixels = width*height;
    for( uint32_t ii=0; ii<numPixels; ++ii )
    {
        if( lastFrame[0] == nextFrame[0] &&
            lastFrame[1] == nextFrame[1] &&
            lastFrame[2] == nextFrame[2] )
        {
            // nextFrame[0] = nextFrame[1] = nextFrame[2] = nextFrame[3] = 0;
            nextFrame[3] = 0;
        }
        else
        {
            nextFrame[3] = 255;
        }
        
        lastFrame += 4;
        nextFrame += 4;
    }
}

void DitherImage( uint8_t* lastFrame, uint8_t* nextFrame, uint32_t width, uint32_t height, Palette& pal )
{
    uint32_t numPixels = width*height;
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
            uint8_t* lastPix = lastFrame? lastFrame + 4*(yy*width+xx) : NULL;
            
            uint8_t rr = nextPix[0] >> 8;
            uint8_t gg = nextPix[1] >> 8;
            uint8_t bb = nextPix[2] >> 8;
            
            if( lastFrame &&
                lastPix[0] == rr &&
                lastPix[1] == gg &&
                lastPix[2] == bb )
            {
                nextPix[0] = rr;
                nextPix[1] = gg;
                nextPix[2] = bb;
                nextPix[3] = 1;
                continue;
            }
            
            int32_t r_err = lastFrame? ((int32_t)nextPix[0]) - (((int32_t)lastPix[0]) << 8) : 1000000;
            int32_t g_err = lastFrame? ((int32_t)nextPix[1]) - (((int32_t)lastPix[1]) << 8) : 1000000;
            int32_t b_err = lastFrame? ((int32_t)nextPix[2]) - (((int32_t)lastPix[2]) << 8) : 1000000;
            
            int32_t bestDiff = abs(r_err)+abs(g_err)+abs(b_err);
            int32_t bestInd = 1;
            
            for( uint32_t jj=0; jj<256; ++jj )
            {
                int32_t r_ierr = ((int32_t)nextPix[0]) - (((int32_t)pal.r[jj]) << 8);
                int32_t g_ierr = ((int32_t)nextPix[1]) - (((int32_t)pal.g[jj]) << 8);
                int32_t b_ierr = ((int32_t)nextPix[2]) - (((int32_t)pal.b[jj]) << 8);
                
                int32_t diff = abs(r_ierr)+abs(g_ierr)+abs(b_ierr);
                if( diff < bestDiff )
                {
                    bestDiff = diff;
                    bestInd = jj;
                    r_err = r_ierr;
                    g_err = g_ierr;
                    b_err = b_ierr;
                }
            }
            
            nextPix[0] = pal.r[bestInd];
            nextPix[1] = pal.g[bestInd];
            nextPix[2] = pal.b[bestInd];
            
            if( lastFrame &&
               lastPix[0] == nextPix[0] &&
               lastPix[1] == nextPix[1] &&
               lastPix[2] == nextPix[2] )
            {
                nextPix[3] = 1;
            }
            else
            {
                nextPix[3] = bestInd;
            }
            
            int quantloc_7 = (yy*width+xx+1);
            int quantloc_3 = (yy*width+width+xx-1);
            int quantloc_5 = (yy*width+width+xx);
            int quantloc_1 = (yy*width+width+xx+1);
            
            if(quantloc_7 < numPixels)
            {
                int32_t* pix7 = quantPixels+4*quantloc_7;
                pix7[0] += r_err * 7 / 16;
                pix7[1] += g_err * 7 / 16;
                pix7[2] += b_err * 7 / 16;
            }
            
            if(quantloc_3 < numPixels)
            {
                int32_t* pix3 = quantPixels+4*quantloc_3;
                pix3[0] += r_err * 3 / 16;
                pix3[1] += g_err * 3 / 16;
                pix3[2] += b_err * 3 / 16;
            }
            
            if(quantloc_5 < numPixels)
            {
                int32_t* pix5 = quantPixels+4*quantloc_5;
                pix5[0] += r_err * 5 / 16;
                pix5[1] += g_err * 5 / 16;
                pix5[2] += b_err * 5 / 16;
            }
            
            if(quantloc_1 < numPixels)
            {
                int32_t* pix1 = quantPixels+4*quantloc_1;
                pix1[0] += r_err / 16;
                pix1[1] += g_err / 16;
                pix1[2] += b_err / 16;
            }
        }
    }
    
    for( uint32_t ii=0; ii<numPixels*4; ++ii )
    {
        nextFrame[ii] = quantPixels[ii];
        if(lastFrame)
            lastFrame[ii] = quantPixels[ii];
    }
    
    free(quantPixels);
}

void ThresholdImage( uint8_t* lastFrame, uint8_t* nextFrame, uint32_t width, uint32_t height, Palette& pal )
{
    uint32_t numPixels = width*height;
    for( uint32_t ii=0; ii<numPixels; ++ii )
    {
        int32_t bestDiff = 1000000;
        int32_t bestInd = 0;
        
        if(lastFrame)
        {
            int32_t r_err = ((int32_t)lastFrame[0]) - ((int32_t)nextFrame[0]);
            int32_t g_err = ((int32_t)lastFrame[1]) - ((int32_t)nextFrame[1]);
            int32_t b_err = ((int32_t)lastFrame[2]) - ((int32_t)nextFrame[2]);
            
            bestDiff = abs(r_err)+abs(g_err)+abs(b_err);
            bestInd = 1;
            
            lastFrame += 4;
        }
        
        for( uint32_t jj=0; jj<256; ++jj )
        {
            int32_t r_err = ((int32_t)pal.r[jj]) - ((int32_t)nextFrame[0]);
            int32_t g_err = ((int32_t)pal.g[jj]) - ((int32_t)nextFrame[1]);
            int32_t b_err = ((int32_t)pal.b[jj]) - ((int32_t)nextFrame[2]);
            
            int32_t diff = abs(r_err)+abs(g_err)+abs(b_err);
            if( diff < bestDiff )
            {
                bestDiff = diff;
                bestInd = jj;
            }
        }
        
        nextFrame[3] = bestInd;
        nextFrame += 4;
    }
}

void CopyToLastFrame( uint8_t* lastFrame, uint8_t* nextFrame, uint32_t width, uint32_t height, Palette& pal )
{
    uint32_t numPixels = width*height;
    for( uint32_t ii=0; ii<numPixels; ++ii )
    {
        int p = nextFrame[3];
        if(p != 1)
        {
            lastFrame[0] = pal.r[p];
            lastFrame[1] = pal.g[p];
            lastFrame[2] = pal.b[p];
        }
        
        lastFrame += 4;
        nextFrame += 4;
    }
}

struct GifWriter
{
    FILE* f;
    Palette pal;
    uint8_t* oldImage;
};

void BeginGif( GifWriter* writer, const char* filename, uint8_t* image, uint32_t width, uint32_t height, uint32_t delay )
{
    writer->f = fopen(filename, "wb");
    ALWAYS_ASSERT(writer->f);
    writer->oldImage = (uint8_t*)malloc(width*height*4);
    
    MakePalette(image, width, height, &writer->pal);
    
    DitherImage(NULL, image, width, height, writer->pal);
    //ThresholdImage(NULL, image, width, height, writer->pal);
    CopyToLastFrame(writer->oldImage, image, width, height, writer->pal);
    
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
    WritePalette(&writer->pal, writer->f);
    
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
    
    WriteLzwImage(writer->f, image, writer->oldImage, 0, 0, width, height, delay, writer->pal);
}

void ContinueGif( GifWriter* writer, uint8_t* image, uint32_t width, uint32_t height, uint32_t delay )
{
    MakePalette(image, width, height, &writer->pal);
    
    DitherImage(writer->oldImage, image, width, height, writer->pal);
    //ThresholdImage(writer->oldImage, image, width, height, writer->pal);
    CopyToLastFrame(writer->oldImage, image, width, height, writer->pal);
    
    WriteLzwImage(writer->f, image, writer->oldImage, 0, 0, width, height, delay, writer->pal);
}

void EndGif( GifWriter* writer )
{
    fputc(0x3b, writer->f); // end of file
    fclose(writer->f);
    free(writer->oldImage);
    
    writer->f = NULL;
    writer->oldImage = NULL;
}

#endif

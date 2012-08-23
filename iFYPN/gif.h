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

typedef std::map<U32, std::pair<U32, U32> > Palette;

U8 Palettize( U8* r, U8* g, U8* b, Palette& pal )
{
    U32 orgb = (*r << 16) | (*g << 8) | *b;
    U32 prgb = pal[orgb].second;
    *r = (prgb >> 16) & 0xff;
    *g = (prgb >> 8) & 0xff;
    *b = (prgb) & 0xff;
    return (prgb >> 24) & 0xff;
}

void AssignColors( Palette* pPal, U32* rPal, U32* gPal, U32* bPal )
{
    for (Palette::iterator it = pPal->begin(); it != pPal->end(); ++it)
    {
        U32 rgb = it->first;
        U32 r = rgb >> 16;
        U32 g = (rgb >> 8) & 0xff;
        U32 b = rgb & 0xff;
        
        U32 bestDiff = 10000;
        U32 bestInd = 256;
        for( U32F ii=0; ii<256; ++ii )
        {
            U32 diff = abs((I32)r-(I32)rPal[ii]) + abs((I32)g-(I32)gPal[ii]) + abs((I32)b-(I32)bPal[ii]);
            if( diff < bestDiff )
            {
                bestDiff = diff;
                bestInd = ii;
            }
        }
        
        U32 prgb = (bestInd << 24) | (rPal[bestInd] << 16) | (gPal[bestInd] << 8) | (bPal[bestInd]);
        it->second.second = prgb;
    }
}

void KMeansIterate( U8* image, U32 width, U32 height, Palette* pPal )
{
    U32 rAvg[256];
    U32 gAvg[256];
    U32 bAvg[256];
    
    U32 count[256];
    
    bzero(rAvg, sizeof(count));
    bzero(gAvg, sizeof(count));
    bzero(bAvg, sizeof(count));
    bzero(count, sizeof(count));
    
    for (Palette::iterator it = pPal->begin(); it != pPal->end(); ++it)
    {
        U32 rgb = it->first;
        U32 r = rgb >> 16;
        U32 g = (rgb >> 8) & 0xff;
        U32 b = rgb & 0xff;
        
        U32 palCount = it->second.first;
        U32 p = it->second.second >> 24;
        
        rAvg[p] += r * palCount;
        gAvg[p] += g * palCount;
        bAvg[p] += b * palCount;
        count[p] += palCount;
    }
    
    for( U32F ii=2; ii<255; ++ii )
    {
        if( count[ii] > 0 )
        {
            U32 bias = count[ii]/2;
            
            rAvg[ii] += bias;
            gAvg[ii] += bias;
            bAvg[ii] += bias;
            
            rAvg[ii] /= count[ii];
            gAvg[ii] /= count[ii];
            bAvg[ii] /= count[ii];
        }
    }
    
    rAvg[0] = 0;
    gAvg[0] = 0;
    bAvg[0] = 0;
    
    rAvg[1] = 0;
    gAvg[1] = 0;
    bAvg[1] = 0;
    
    rAvg[255] = 255;
    gAvg[255] = 255;
    bAvg[255] = 255;
    
    std::vector< std::pair<U64, U32> > paletteDists;
    for (Palette::iterator it = pPal->begin(); it != pPal->end(); ++it)
    {
        U32 orgb = it->first;
        I32 oor = orgb >> 16;
        I32 og = (orgb >> 8) & 0xff;
        I32 ob = orgb & 0xff;
        
        U32 prgb = it->second.second;
        I32 pr = prgb >> 16;
        I32 pg = (prgb >> 8) & 0xff;
        I32 pb = prgb & 0xff;
        
        U32 pcount = it->second.first;
        U64 dist = (abs(oor-pr) + abs(og-pg) + abs(ob-pb)) * pcount;
        
        paletteDists.push_back(std::pair<U64, U32>(dist, orgb));
    }
    std::sort(paletteDists.begin(), paletteDists.end());
    
    for( U32F ii=2; ii<255; ++ii )
    {
        if( count[ii] == 0 && paletteDists.size() > 0 )
        {
            U64 dist = paletteDists.back().first;
            U32 col = paletteDists.back().second;
            paletteDists.pop_back();
            
            rAvg[ii] = col >> 16;
            gAvg[ii] = (col >> 8) & 0xff;
            bAvg[ii] = col & 0xff;
        }
    }
    
    AssignColors(pPal, rAvg, gAvg, bAvg);
}

void MakePalette( U8* image, U32 width, U32 height, Palette* pPal )
{
    U32 numPixels = width*height;
    for( U32F ii=0; ii<numPixels; ++ii )
    {
        U8 r = image[ii*4];
        U8 g = image[ii*4+1];
        U8 b = image[ii*4+2];
        
        U32 rgb = (r << 16) | (g << 8) | b;
        
        if(pPal->find(rgb) != pPal->end())
        {
            (*pPal)[rgb].first++;
        }
        else
        {
            (*pPal)[rgb] = std::pair<U32,U32>(0,0);
        }
    }
    
    U32 rBas[256];
    U32 gBas[256];
    U32 bBas[256];
    
    rBas[0] = 0;
    gBas[0] = 0;
    bBas[0] = 0;
    
    rBas[1] = 0;
    gBas[1] = 0;
    bBas[1] = 0;
    
    rBas[255] = 255;
    gBas[255] = 255;
    bBas[255] = 255;
    
    for( U32F ii=2; ii<255; ++ii )
    {
        U32 x = rand() % width;
        U32 y = rand() % height;
        U8 r = image[(y*width+x)*4];
        U8 g = image[(y*width+x)*4+1];
        U8 b = image[(y*width+x)*4+2];
        rBas[ii] = r;
        gBas[ii] = g;
        bBas[ii] = b;
    }
    
    AssignColors(pPal, rBas, gBas, bBas);
    
    for(U32F ii=0; ii<5; ++ii)
    {
        KMeansIterate(image, width, height, pPal);
    }
}

struct BitStatus
{
    U8 bitIndex;
    U8 byte;
    
    U32 chunkIndex;
    U8 chunk[256];
};

void WriteBit( BitStatus& stat, U32 bit )
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

void WriteCode( FILE* f, BitStatus& stat, U32 code, U32 length )
{
    //U32 bitsRemaining = (255 - stat.chunkIndex)*8 + 8-stat.bitIndex;
    
    //if( bitsRemaining < length )
    //{
    //    WriteImageChunk(f, stat);
    //}
    
    for( U32F ii=0; ii<length; ++ii )
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
    U32 m_code;
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
    
    for( U32F ii=0; ii<256; ++ii )
        DeleteLzwTree(node->m_children[ii]);
    
    free(node);
}

void WritePalette( const Palette* pPal, FILE* f )
{
    U32 cols[256];
    
    for (Palette::const_iterator it = pPal->begin(); it != pPal->end(); ++it)
    {
        U32 prgb = it->second.second;
        U32 p = (prgb >> 24) & 0xff;
        
        cols[p] = prgb;
    }
    
    fputc(0, f);  // first and second colors: always black
    fputc(0, f);
    fputc(0, f);
    
    fputc(0, f);
    fputc(0, f);
    fputc(0, f);
    
    for(U32F ii=2; ii<255; ++ii)
    {
        U32 r = (cols[ii] >> 16) & 0xff;
        U32 g = (cols[ii] >> 8) & 0xff;
        U32 b = (cols[ii]) & 0xff;
        
        fputc(r, f);
        fputc(g, f);
        fputc(b, f);
    }
    fputc(255, f);  // last color: always white
    fputc(255, f);
    fputc(255, f);
}

void WriteLzwImage(FILE* f, U8* image, U8* oldImage, U32 left, U32 top,  U32 width, U32 height, U32 delay, Palette& pal)
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
    
    //fputc(0, f); // no local color table
    fputc(0x87, f); // local color table present, 256 entries
    
    WritePalette(&pal, f);
    
    fputc(8, f); // min code size 8 bits
    
    lzwNode* codetree = InitLzwTree();
    lzwNode* curNode = codetree;
    U32 codeSize = 9;
    U32 maxCode = 257;
    
    BitStatus stat;
    stat.byte = 0;
    stat.bitIndex = 0;
    stat.chunkIndex = 0;
    
    WriteCode(f, stat, 256, codeSize);
    
    for(I32F yy=height-1; yy>=0; --yy)
    {
        for(I32F xx=0; xx<width; ++xx)
        {
            U8 r = image[(yy*width+xx)*4];
            U8 g = image[(yy*width+xx)*4+1];
            U8 b = image[(yy*width+xx)*4+2];
            U8 a = image[(yy*width+xx)*4+4];
            U8 nextValue = Palettize(&r, &g, &b, pal);
            
            if( a == 0 )
            {
                nextValue = 1;
            }
            else
            {
                oldImage[(yy*width+xx)*4] = r;
                oldImage[(yy*width+xx)*4+1] = g;
                oldImage[(yy*width+xx)*4+2] = b;
            }
            
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
                    if( codeSize == 12 )
                    {
                        WriteCode(f, stat, 256, codeSize); // clear tree
                        
                        DeleteLzwTree(codetree);
                        codetree = InitLzwTree();
                        curNode = codetree;
                        codeSize = 9;
                        maxCode = 257;
                    }
                }
                
                curNode = codetree->m_children[nextValue];
            }
        }
    }
    
    WriteCode(f, stat, 257, codeSize);
    while( stat.bitIndex ) WriteBit(stat, 0);
    if( stat.chunkIndex ) WriteImageChunk(f, stat);
    DeleteLzwTree(codetree);
    fputc(0, f); // image block terminator
}

void SetTransparency( U8* oldImg, U8* newImg, U32 width, U32 height )
{
    U32 numPixels = width*height;
    for( U32F ii=0; ii<numPixels; ++ii )
    {
        if( oldImg[0] == newImg[0] &&
            oldImg[1] == newImg[1] &&
            oldImg[2] == newImg[2] )
            newImg[3] = 0;
        else
            newImg[3] = 255;
        
        oldImg += 4;
        newImg += 4;
    }
}

void WriteOverTransparency( U8* img, U32 width, U32 height )
{
    U32 numPixels = width*height;
    for( U32F ii=0; ii<numPixels; ++ii )
    {
        img[3] = 255;
        img += 4;
    }
}

struct GifWriter
{
    FILE* f;
    U8* oldImage;
};

void BeginGif( GifWriter* writer, U8* image, U32 width, U32 height, U32 delay )
{
    writer->f = fopen("/Users/ctangora/cmt-test.gif", "wb");
    ALWAYS_ASSERT(writer->f);
    
    Palette pal;
    MakePalette(image, width, height, &pal);
    WriteOverTransparency(image, width, height);
    
    writer->oldImage = (U8*)malloc(width*height*4);
    memcpy(writer->oldImage, image, width*height*4);
    
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
    WritePalette(&pal, writer->f);
    
    if( delay != 0 )
    {
        // animation header
        fputc(0x21, writer->f); // extension
        fputc(0xff, writer->f); // application specific
        fputc(11, writer->f); // length 11
        fputs("NETSCAPE2.0", writer->f); // yes, really
        fputc(3, writer->f); // 3 bytes of NETSCAPE2.0 data
        
        fputc(0, writer->f); // loop infinitely (byte 0)
        fputc(0, writer->f); // loop infinitely (byte 1)
        
        fputc(0, writer->f); // block terminator
    }
    
    WriteLzwImage(writer->f, image, writer->oldImage, 0, 0, width, height, delay, pal);
}

void ContinueGif( GifWriter* writer, U8* image, U32 width, U32 height, U32 delay )
{
    Palette pal;
    MakePalette(image, width, height, &pal);
    SetTransparency(writer->oldImage, image, width, height);
    
    WriteLzwImage(writer->f, image, writer->oldImage, 0, 0, width, height, delay, pal);
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

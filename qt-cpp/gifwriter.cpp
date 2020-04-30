#include "gifwriter.h"
#include <QBuffer>
#include <QImage>
#include <QSize>

namespace
{

using RGB = MyLib::GifWriter::RGB;

RGB operator+(const RGB& a, const RGB& b)
{
  RGB result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

RGB operator-(const RGB& a, const RGB& b)
{
  RGB result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

RGB& operator+=(RGB& a, const RGB& b)
{
  a = a + b;
  return a;
}

RGB operator/(const RGB& rgb, double f)
{
  RGB result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = rgb[i] / f;
  }
  return result;
}

RGB operator*(const RGB& rgb, double f)
{
  RGB result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = rgb[i] * f;
  }
  return result;
}

RGB min(const RGB& a, const RGB& b)
{
  RGB result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = std::min(a[i], b[i]);
  }
  return result;
}

RGB max(const RGB& a, const RGB& b)
{
  RGB result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = std::max(a[i], b[i]);
  }
  return result;
}

template<typename T> RGB rgb_from_ptr(const T* ptr)
{
  return { ptr[0], ptr[1], ptr[2] };
}

template<typename T> void rgb_to_ptr(T* ptr, const RGB& rgb)
{
  for (std::size_t i = 0; i < 3; ++i) {
    ptr[i] = rgb[i];
  }
}

RGB operator-(const RGB& rgb)
{
  return { -rgb[0], -rgb[1], -rgb[2] };
}

void swap_pixels(uint8_t* image, int pixA, int pixB)
{
  for (std::size_t i = 0; i < 4; ++i) {
    std::swap(image[4*pixA+i], image[4*pixB+i]);
  }
}

// just the partition operation from quicksort
int partition(uint8_t* image, const int left, const int right, const int elt, int pivotIndex)
{
    const int pivotValue = image[(pivotIndex)*4+elt];
    swap_pixels(image, pivotIndex, right-1);
    int storeIndex = left;
    bool split = 0;
    for (int ii = left; ii < right-1; ++ii) {
      int arrayVal = image[ii*4+elt];
      if (arrayVal < pivotValue) {
        swap_pixels(image, ii, storeIndex);
        ++storeIndex;
      } else if (arrayVal == pivotValue) {
        if (split) {
          swap_pixels(image, ii, storeIndex);
          ++storeIndex;
        }
        split = !split;
      }
    }
    swap_pixels(image, storeIndex, right-1);
    return storeIndex;
}

// Perform an incomplete sort, finding all elements above and below the desired median
void partition_by_median(uint8_t* image, int left, int right, int com, int neededCenter)
{
  if (left < right-1) {
    int pivotIndex = left + (right-left)/2;

    pivotIndex = partition(image, left, right, com, pivotIndex);

    // Only "sort" the section of the array that contains the median
    if (pivotIndex > neededCenter) {
      partition_by_median(image, left, pivotIndex, com, neededCenter);
    }

    if (pivotIndex < neededCenter) {
      partition_by_median(image, pivotIndex+1, right, com, neededCenter);
    }
  }
}

// Finds all pixels that have changed from the previous image and
// moves them to the fromt of th buffer.
// This allows us to build a palette optimized for the colors of the
// changed pixels only.
int pick_changed_pixels(const uint8_t* last_frame, uint8_t* frame, int num_pixels)
{
  int numChanged = 0;
  uint8_t* writeIter = frame;

  for (int ii = 0; ii < num_pixels; ++ii) {
    if (last_frame[0] != frame[0] || last_frame[1] != frame[1] || last_frame[2] != frame[2]) {
      writeIter[0] = frame[0];
      writeIter[1] = frame[1];
      writeIter[2] = frame[2];
      ++numChanged;
      writeIter += 4;
    }
    last_frame += 4;
    frame += 4;
  }
  return numChanged;
}

}  // namespace


namespace MyLib
{

GifWriter::GifWriter(const QSize& size, int delay, QBuffer& buffer)
  : m_size(size), m_delay(delay), m_buffer(buffer)
  , m_old_image(size.width(), size.height(), image_format)
{
  m_old_image.fill(0);
  puts("GIF89a");

  // screen descriptor
  {
    const uint32_t width = m_size.width();
    const uint32_t height = m_size.height();
    putc(width & 0xff);
    putc((width >> 8) & 0xff);
    putc(height & 0xff);
    putc((height >> 8) & 0xff);
  }

  putc(0xf0);
  putc(0);     // background color
  putc(0);     // pixels are square (we need to specify this because it's 1989)

  // now the "global" palette (really just a dummy palette)
  put(RGB{0, 0, 0});  // color 0: black
  put(RGB{0, 0, 0});  // color 1: also black

  if (delay != 0) {
    // animation header
    putc(0x21); // extension
    putc(0xff); // application specific
    putc(11); // length 11
    puts("NETSCAPE2.0"); // yes, really
    putc(3); // 3 bytes of NETSCAPE2.0 data

    putc(1); // JUST BECAUSE
    putc(0); // loop infinitely (byte 0)
    putc(0); // loop infinitely (byte 1)

    putc(0); // block terminator
  }
}

void GifWriter::write_frame(QImage image)
{
  image = image.convertToFormat(image_format).scaled(m_size);

  uint8_t* old_image_bits = m_first_frame ? nullptr : m_old_image.bits();
  const GifPalette palette(dither || m_first_frame ? QImage() : m_old_image, image, dither);

  if (dither) {
    dither_image(old_image_bits, image.bits(), m_old_image.bits(), palette);
  } else {
    threshold_image(old_image_bits, image.bits(), m_old_image.bits(), palette);
  }

  write_lzw_image(palette);
  m_first_frame = false;
}

QByteArray GifWriter::encode(const std::vector<QImage>& frames, int delay, bool dither)
{
  QByteArray data;
  if (frames.size() > 0) {
    const QSize size = frames.front().size();
    QBuffer buffer(&data);
    buffer.open(QIODevice::WriteOnly);
    assert(buffer.isWritable());
    GifWriter writer(frames.front().size(), delay, buffer);
    writer.dither = dither;
    for (const QImage& image : frames) {
      assert(image.size() == size);
      writer.write_frame(image);
    }
  }
  return data;
}

void GifWriter::puts(const char* str)
{
  puts(str, strlen(str));
}

void GifWriter::puts(const char* str, std::size_t n)
{
  const std::size_t m = m_buffer.write(str, n);
  assert(n == m);
  Q_UNUSED(m)
}

void GifWriter::putc(const char c)
{
  m_buffer.putChar(c);
}

void GifWriter::put(const GifWriter::RGB& rgb)
{
  for (std::size_t i = 0; i < rgb.size(); ++i) {
    putc(rgb[i]);
  }
}

void GifWriter::dither_image(const uint8_t* last_frame, const uint8_t* next_frame,
                             uint8_t* out_frame, const GifWriter::GifPalette& palette)
{
  const int num_pixels = m_size.width() * m_size.height();

  // quantPixels initially holds color*256 for all pixels
  // The extra 8 bits of precision allow for sub-single-color error values
  // to be propagated
  std::vector<int32_t> quant_pixels(num_pixels * 4, 0);

  for (int ii = 0; ii < num_pixels * 4; ++ii) {
    uint8_t pix = next_frame[ii];
    int32_t pix16 = int32_t(pix) * 256;
    quant_pixels[ii] = pix16;
  }

  const int width = m_size.width();
  const int height = m_size.height();

  for (int yy = 0; yy < height; ++yy) {
    for (int xx = 0; xx < width; ++xx) {
      const std::size_t offset = 4 * (yy * width + xx);
      int32_t* next_pix = quant_pixels.data() + offset;

      // Compute the colors we want (rounding to nearest)
      RGB rrggbb = (rgb_from_ptr(next_pix) + RGB {127, 127, 127}) / 256;

      // if it happens that we want the color from last frame, then just write out
      // a transparent pixel
      if (!m_old_image.isNull()) {
        if (last_frame != nullptr && rgb_from_ptr(&last_frame[offset]) == rrggbb) {
          for (std::size_t i = 0; i < rrggbb.size(); ++i) {
            next_pix[i] = rrggbb[i];
          }
          next_pix[3] = k_gif_trans_index;
          continue;
        }
      }

      int32_t best_diff = 1000000;
      int32_t best_index = k_gif_trans_index;

      // Search the palete
      get_closest_palette_color(palette, rrggbb, best_index, best_diff);

      // Write the result to the temp buffer
      const RGB err = rgb_from_ptr(next_pix) - (palette.rgb[best_index] * 256);

      for (std::size_t i = 0; i < err.size(); ++i) {
        next_pix[i] = palette.rgb[best_index][i];
      }
      next_pix[3] = best_index;

      int32_t* data = quant_pixels.data();
      auto bm = [data, num_pixels](int quantloc, RGB err, int denom) {
        if (quantloc < num_pixels) {
          int32_t* pix = data + 4 * quantloc;
          rgb_to_ptr(pix, rgb_from_ptr(pix) + max(-rgb_from_ptr(pix), err * denom / 16));
        }
      };

//      // Propagate the error to the four adjacent locations
//      // that we haven't touched yet
      bm(yy * width         + xx + 1, err, 7);
      bm(yy * width + width + xx - 1, err, 3);
      bm(yy * width + width + xx,     err, 5);
      bm(yy * width + width + xx + 1, err, 1);
    }
  }

  // Copy the palettized result to the output buffer
  for (int ii = 0; ii < num_pixels * 4; ++ii) {
    out_frame[ii] = quant_pixels[ii];
  }
}

void GifWriter::threshold_image(const uint8_t* last_frame, const uint8_t* next_frame,
                                uint8_t* out_frame, const GifWriter::GifPalette& palette)
{
  uint32_t num_pixels = m_size.width() * m_size.height();

  for (uint32_t ii = 0; ii < num_pixels; ++ii) {
    // if a previous color is available, and it matches the current color,
    // set the pixel to transparent

    if (last_frame && rgb_from_ptr(last_frame) == rgb_from_ptr(next_frame)) {
      rgb_to_ptr(out_frame, rgb_from_ptr(last_frame));
      out_frame[3] = k_gif_trans_index;
    } else {
      // palettize the pixel
      int32_t best_diff = 1000000;
      int32_t best_index = 1;
      get_closest_palette_color(palette,
                                { next_frame[0], next_frame[1], next_frame[2] },
                                best_index, best_diff);

      // Write the resulting color to the output buffer
      rgb_to_ptr(out_frame, palette.rgb[best_index]);
      out_frame[3] = best_index;
    }

    if (last_frame) {
      last_frame += 4;
    }

    out_frame += 4;
    next_frame += 4;
  }

}

void GifWriter::write_lzw_image(const GifWriter::GifPalette& palette)
{
  assert(!m_old_image.isNull());
  const uint8_t* image_bits = m_old_image.bits();

  // graphics control extension
  putc(0x21);
  putc(0xf9);
  putc(0x04);
  putc(0x05); // leave prev frame in place, this frame has transparency
  putc(m_delay & 0xff);
  putc((m_delay >> 8) & 0xff);
  putc(k_gif_trans_index); // transparent color index
  putc(0);

  putc(0x2c); // image descriptor block

  const uint32_t left = 0;
  const uint32_t top = 0;
  const uint32_t width = m_size.width();
  const uint32_t height = m_size.height();

  putc(left & 0xff);           // corner of image in canvas space
  putc((left >> 8) & 0xff);
  putc(top & 0xff);
  putc((top >> 8) & 0xff);

  putc(width & 0xff);          // width and height of image
  putc((width >> 8) & 0xff);
  putc(height & 0xff);
  putc((height >> 8) & 0xff);

  //fputc(0); // no local color table, no transparency
  //fputc(0x80); // no local color table, but transparency

  putc(0x80 + palette.bit_depth-1); // local color table present, 2 ^ bitDepth entries
  write_palette(palette);

  const int minCodeSize = palette.bit_depth;
  const uint32_t clearCode = 1 << palette.bit_depth;

  putc(minCodeSize); // min code size 8 bits

  std::vector<GifLzwNode> codetree(4096);
  memset(codetree.data(), 0, sizeof(GifLzwNode)*4096);

  int32_t curCode = -1;
  uint32_t codeSize = (uint32_t)minCodeSize + 1;
  uint32_t maxCode = clearCode+1;

  GifBitStatus stat;
  stat.byte = 0;
  stat.bit_index = 0;
  stat.chunk_index = 0;

  write_code(stat, clearCode, codeSize);  // start with a fresh LZW dictionary

  for (uint32_t yy = 0; yy < height; ++yy) {
    for (uint32_t xx = 0; xx < width; ++xx) {
#ifdef GIF_FLIP_VERT
      // bottom-left origin image (such as an OpenGL capture)
      uint8_t nextValue = image_bits[((height-1-yy)*width+xx)*4+3];
#else
      // top-left origin
      uint8_t nextValue = image_bits[(yy*width+xx)*4+3];
#endif

      // "loser mode" - no compression, every single code is followed immediately by a clear
      //WriteCode( f, stat, nextValue, codeSize );
      //WriteCode( f, stat, 256, codeSize );

      if (curCode < 0) {
        // first value in a new run
        curCode = nextValue;
      } else if (codetree[curCode].m_next[nextValue]) {
        // current run already in the dictionary
        curCode = codetree[curCode].m_next[nextValue];
      } else {
        // finish the current run, write a code
        write_code(stat, (uint32_t)curCode, codeSize);

        // insert the new run into the dictionary
        codetree[curCode].m_next[nextValue] = (uint16_t)++maxCode;

        if (maxCode >= (1ul << codeSize)) {
          // dictionary entry count has broken a size barrier,
          // we need more bits for codes
          codeSize++;
        }

        if (maxCode == 4095) {
          // the dictionary is full, clear it out and begin anew
          write_code(stat, clearCode, codeSize); // clear tree

          memset(codetree.data(), 0, sizeof(GifLzwNode)*4096);
          codeSize = (uint32_t)(minCodeSize + 1);
          maxCode = clearCode+1;
        }

        curCode = nextValue;
      }
    }
  }

  // compression footer
  write_code(stat, curCode, codeSize);
  write_code(stat, clearCode, codeSize);
  write_code(stat, clearCode + 1, minCodeSize + 1);

  // write out the last partial chunk
  while (stat.bit_index) {
    write_bit(stat, 0);
  }

  if (stat.chunk_index) {
    write_chunk(stat);
  }

  putc(0); // image block terminator
}

GifWriter::GifPalette::GifPalette(const QImage& last_frame, const QImage& next_frame, bool dither)
  : bit_depth(8), rgb(256, {0, 0, 0}), tree_split_elt(256, 0), tree_split(256, 0), m_dither(dither)
{

  const int width = next_frame.width();
  const int height = next_frame.height();
  assert(last_frame.isNull() || last_frame.size() == next_frame.size());

  // SplitPalette is destructive (it sorts the pixels by color) so
  // we must create a copy of the image for it to destroy
  QImage destroyable_image = next_frame.copy();

  int num_pixels = width * height;
  if (!last_frame.isNull()) {
    num_pixels = pick_changed_pixels(last_frame.bits(), destroyable_image.bits(), num_pixels);
  }

  const int last_elt = 1 << bit_depth;
  const int split_elt = last_elt/2;
  const int split_dist = split_elt/2;
  split(destroyable_image.bits(), num_pixels, 1, last_elt, split_elt, split_dist, 1);

  // add the bottom node for the transparency index
  tree_split[1 << (bit_depth-1)] = 0;
  tree_split_elt[1 << (bit_depth-1)] = 0;

  rgb[0] = {0, 0, 0};
}

void GifWriter::write_code(GifBitStatus& stat, uint32_t code, uint32_t length)
{
  for (uint32_t ii = 0; ii < length; ++ii) {
    write_bit(stat, code);
    code = code >> 1;

    if (stat.chunk_index == 255) {
      write_chunk(stat);
    }
  }
}

void GifWriter::write_chunk(GifBitStatus& stat)
{
  putc(stat.chunk_index);
  puts(reinterpret_cast<const char*>(stat.chunk.data()), 1 * stat.chunk_index);

  stat.bit_index = 0;
  stat.byte = 0;
  stat.chunk_index = 0;
}

void GifWriter::write_bit(GifBitStatus& stat, uint32_t bit)
{
  bit = bit & 1;
  bit = bit << stat.bit_index;
  stat.byte |= bit;

  ++stat.bit_index;
  if (stat.bit_index > 7) {
    // move the newly-finished byte to the chunk buffer
    stat.chunk[stat.chunk_index++] = stat.byte;
    // and start a new byte
    stat.bit_index = 0;
    stat.byte = 0;
  }
}


// write a 256-color (8-bit) image palette to the file
void GifWriter::write_palette(const GifWriter::GifPalette& palette)
{
  put(RGB {0, 0, 0});
  for (int ii = 1; ii < (1 << palette.bit_depth); ++ii) {
    put(palette.rgb[ii]);
  }
}

// Writes the EOF code, closes the file handle, and frees temp memory used by a GIF.
// Many if not most viewers will still display a GIF properly if the EOF code is missing,
// but it's still a good idea to write it out.
GifWriter::~GifWriter()
{
  putc(0x3b); // end of file
}

// walks the k-d tree to pick the palette entry for a desired color.
// Takes as in/out parameters the current best color and its error -
// only changes them if it finds a better color in its subtree.
// this is the major hotspot in the code at the moment.
void GifWriter::get_closest_palette_color(const GifPalette& palette, const RGB& rgb,
                                          int& best_index, int& best_diff, int treeRoot)
{
  // base case, reached the bottom of the tree
  if (treeRoot > (1 << palette.bit_depth)-1) {
    int ind = treeRoot - (1 << palette.bit_depth);
    if (ind == k_gif_trans_index) {
      return;
    }

    // check whether this color is better than the current winner

    int diff = 0;
    for (std::size_t i = 0; i < 3; ++i) {
      diff += std::abs(rgb[i] - palette.rgb[ind][i]);
    }

    if (diff < best_diff) {
      best_index = ind;
      best_diff = diff;
    }

    return;
  }

  // take the appropriate color (r, g, or b) for this node of the k-d tree
  int splitComp = rgb[palette.tree_split_elt[treeRoot]];

  int splitPos = palette.tree_split[treeRoot];
  if (splitPos > splitComp) {
      // check the left subtree
      get_closest_palette_color(palette, rgb, best_index, best_diff, treeRoot*2);
      if (best_diff > splitPos - splitComp) {
        // cannot prove there's not a better value in the right subtree, check that too
        get_closest_palette_color(palette, rgb, best_index, best_diff, treeRoot*2+1);
      }
  } else {
    get_closest_palette_color(palette, rgb, best_index, best_diff, treeRoot*2+1);
    if (best_diff > splitComp - splitPos) {
      get_closest_palette_color(palette, rgb, best_index, best_diff, treeRoot*2);
    }
  }
}

// Builds a palette by creating a balanced k-d tree of all pixels in the image
void GifWriter::GifPalette::split(uint8_t* image, int num_pixels, int first_elt, int last_elt, int split_elt, int split_dist, int treeNode)
{
    if (last_elt <= first_elt || num_pixels == 0) {
      return;
    }

    // base case, bottom of the tree
    if (last_elt == first_elt+1) {
      if (m_dither) {
        // Dithering needs at least one color as dark as anything
        // in the image and at least one brightest color -
        // otherwise it builds up error and produces strange artifacts
        if (first_elt == 1) {
          // special case: the darkest color in the image
          RGB rgb = {255, 255, 255};
          for (int ii = 0; ii < num_pixels; ++ii) {
            rgb = min(rgb, rgb_from_ptr(&image[ii * 4]));
          }
          this->rgb[first_elt] = rgb;
          return;
        }
        if (first_elt == (1 << bit_depth)-1) {
          // special case: the lightest color in the image
          RGB rgb {0, 0, 0};
          for (int ii = 0; ii < num_pixels; ++ii) {
            rgb = max(rgb, rgb_from_ptr(&image[ii*4]));
          }
          this->rgb[first_elt] = rgb;
          return;
        }
      }

      // otherwise, take the average of all colors in this subcube
      RGB rgb = {0, 0, 0};
      for (int ii=0; ii<num_pixels; ++ii) {
        rgb += rgb_from_ptr(&image[ii*4]);
      }

      rgb = (rgb + RGB{ num_pixels, num_pixels, num_pixels } / 2) / num_pixels;  // round to nearest

      this->rgb[first_elt] = rgb;

      return;
    }

    // Find the axis with the largest range
    RGB min_rgb {255, 255, 255};
    RGB max_rgb {0, 0, 0};
    for (int ii = 0; ii < num_pixels; ++ii) {
      RGB rgb = rgb_from_ptr(&image[ii*4]);
      max_rgb = max(rgb, max_rgb);
      min_rgb = min(rgb, min_rgb);
    }

    RGB rgb_range = max_rgb - min_rgb;

    // and split along that axis. (incidentally, this means this isn't a "proper" k-d tree but I don't know what else to call it)
    int splitCom = 1;
    if (rgb_range[2] > rgb_range[1]) {
      splitCom = 2;
    }
    if (rgb_range[0] > rgb_range[2] && rgb_range[0] > rgb_range[1]) {
      splitCom = 0;
    }

    int sub_pixel_a = num_pixels * (split_elt - first_elt) / (last_elt - first_elt);
    int sub_pixel_b = num_pixels-sub_pixel_a;

    partition_by_median(image, 0, num_pixels, splitCom, sub_pixel_a);

    tree_split_elt[treeNode] = splitCom;
    tree_split[treeNode] = image[sub_pixel_a*4+splitCom];

    split(image,               sub_pixel_a, first_elt, split_elt, split_elt-split_dist, split_dist/2, treeNode*2);
    split(image+sub_pixel_a*4, sub_pixel_b, split_elt, last_elt,  split_elt+split_dist, split_dist/2, treeNode*2+1);
}

GifWriter::GifBitStatus::GifBitStatus()
  : chunk(256, 0)
{
}

}  // namespace MyLib

#pragma once

// Based on code by Charlie Tangora (ctangora@gmail.com) (Public domain)

#include <QImage>
#include <QSize>
#include <QDataStream>
#include <QImage>
#include <array>

class QBuffer;

namespace MyLib
{

class GifWriter
{
public:
  using RGB = std::array<int, 3>;
  GifWriter(const QSize& size, int delay, QBuffer& buffer);
  ~GifWriter();
  void write_frame(QImage image);
  static QByteArray encode(const std::vector<QImage>& frames, int delay, bool dither);
  bool dither = true;
  static constexpr QImage::Format image_format = QImage::Format_RGBA8888_Premultiplied;

private:
  const QSize m_size;
  int m_delay = 200;
  QBuffer& m_buffer;
  QImage m_old_image;
  bool m_first_frame = true;
  static constexpr int k_gif_trans_index = 0;

  void puts(const char* str);
  void puts(const char* str, std::size_t n);
  void putc(const char c);
  void put(const RGB& rgb);
  void write_frame(const uint8_t* data);

  // Simple structure to write out the LZW-compressed portion of the image
  // one bit at a time
  struct GifBitStatus
  {
    GifBitStatus();
    uint8_t bit_index; // how many bits in the partial byte written so far
    uint8_t byte;      // current partial byte

    uint32_t chunk_index;

    // bytes are written in here until we have 256 of them, then written to the file
    std::vector<uint8_t> chunk;
  };

  struct GifPalette
  {
    GifPalette(const QImage& last_frame, const QImage& next_frame, bool dither);
    int bit_depth;

    std::vector<RGB> rgb;

    // k-d tree over RGB space, organized in heap fashion
    // i.e. left child of node i is node i*2, right child is node i*2+1
    // nodes 256-511 are implicitly the leaves, containing a color
    std::vector<uint8_t> tree_split_elt;
    std::vector<uint8_t> tree_split;
    void split(uint8_t* image, int numPixels, int firstElt, int lastElt, int splitElt, int splitDist, int treeNode);

  private:
    const bool m_dither;

  };

  // The LZW dictionary is a 256-ary tree constructed as the file is encoded,
  // this is one node
  struct GifLzwNode
  {
    uint16_t m_next[256];
  };

  void dither_image(const uint8_t* last_frame, const uint8_t* next_frame, uint8_t* out_frame, const GifPalette& palette);
  void threshold_image(const uint8_t* last_frame, const uint8_t* next_frame, uint8_t* out_frame, const GifPalette& palette);
  void write_lzw_image(const GifPalette& palette);
  void write_code(GifBitStatus& stat, uint32_t code, uint32_t length);
  void write_chunk(GifBitStatus& stat);
  void write_bit(GifBitStatus& stat, uint32_t bit);
  void write_palette(const GifPalette& palette);
  void get_closest_palette_color(const GifPalette& palette, const RGB& rgb,
                                 int& best_index, int& bestDiff, int treeRoot = 1);
};

}  // namespace MyLib
